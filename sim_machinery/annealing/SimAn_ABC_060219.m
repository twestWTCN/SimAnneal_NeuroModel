function [R,parBank] = SimAn_ABC_060219(R,p,m,parBank)
%%%% SIMULATED ANNEALING for APROXIMATE BAYESIAN COMPUTATION for
%%%% HIGH DIMENSIONAL DYNAMICAL MODELS
% ---- 16/08/18---------------------------
% This annealing function uses approximate Bayesian computation in order to
% estimate the posterior parameter distributions, using a shifting epsilon
% that moves with a cooling schedule.
% Notes
% This version collapses parameter resamples into a single function
% adapted for generic usage
% ic - starting conditions
% u - external input (can be empty
% p - structure of parameters
% m - structure of model specfications
% R - settings for annealing, plotting, integration etc.
%
% This function will output R at each annealing loop (assigned into 'base'
% workspace. The field R.mfit is appended with the Rho (correlation matrix) and nu
% (degrees of freedom) of the estimated copula. Multivariate Gaussian
% is also estimated and specified in R.mfit.Mu and R.mfit.Sigma.
% There are several plotting functions which will track the progress of the
% annealing.
% TO DO:
% Change convergence parameters upon the inference of parameter precision.
% In this way the convergence depends upon the inference and not on the
% outcome (accuracy).
% Timothy West (2018) - UCL CoMPLEX
% / UCL, Wellcome Trust Centre for Human Neuroscience
%%%%%%%%%%%%%%%%%%%%%%
%% Setup for annealing
if nargin<4
    parBank = [];
end
if ~isfield(R.plot,'flag')
    R.plot.flag = 1; % plotting is default behaviour
end
if isempty(m)
    m.m = 1;
end
pOrg = p; % Record prior parameters.

% Set Fixed Annealing Parameters
eps_prior = -4; % prior eps (needed for gradient approximation);
eps_exp = -3.9;
eps_act = eps_prior;
delta_act = 0.05;
% Compute indices of parameters to be optimized
[pInd,pMu,pSig] = parOptInds_110817(R,p,m.m); % in structure form
% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
pMuMap = spm_vec(pMu);
pSigMap = spm_vec(pSig);
R.SimAn.minRank = ceil(size(pIndMap,1)*3); %Ensure rank of sample is large enough to compute copula
% set initial batch of parameters from gaussian priors
if isfield(R,'Mfit')
    rep =  R.SimAn.rep(1);
    par = postDrawCopula(R.Mfit,p,pIndMap,rep);
else
    rep = R.SimAn.rep(1);
    ptmp = spm_vec(p);
    Mfit.Mu = ptmp(pMuMap);
    Mfit.Sigma = ptmp(pSigMap).*R.SimAn.jitter;
    par = postDrawMVN(Mfit,pOrg,pIndMap,pSigMap,rep);
end
itry = 0; cflag = 0;
ii = 1;
%% Main Annealing Loop
while ii <= R.SimAn.searchMax
    %% Batch Loop for Replicates for Generation of Pseudodata
    % This is where the heavy work is done. This is run inside parfor. Any
    % optimization here is prime.
    clear xsims_rep feat_sim_rep
    parfor jj = 1:rep % Replicates for each temperature
        % Get sample Parameters
        pnew = par{jj};
        %% Simulate New Data
        u = innovate_timeseries(R,m);
        u{1} = u{1}.*sqrt(R.IntP.dt);
        [r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData(R,m,u,pnew,0,0);
        
        r2rep{jj} = r2;
        par_rep{jj} = pnew;
        xsims_rep{jj} = xsims_gl; % This takes too much memory: !Modified to store last second only!
        feat_sim_rep{jj} = feat_sim;
        disp(['Iterate ' num2str(jj) ' temperature ' num2str(ii)])
    end % End of batch replicates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARAMETER OPTIMIZATION BEGINS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % icop(1)>itry: Try to form copula with Annealing or Percentile eps
    % icop(2)>itry>icop(1): Find eps to form minimum rank from parbank
    % icop(2)<itry: Try to force
    
    %% Concatanate Batch Results and Decide Acceptance Level Epsilon
    % Retrieve fits
    r2loop = [r2rep{:}];
    % Delete failed simulations
    r2loop(r2loop==1) = -inf;
    r2loop(isnan(r2loop)==1) = -inf;
    r2loop(imag(r2loop)==1) = -inf;
    r2loop(isinf(r2loop)==1) = -inf;
    % Append succesful replicates to bank of params and fits
    %(parameter table, with fits)
    for i = 1:numel(r2loop)
        if ~isinf(r2loop(i))
            parI(:,i) = [full(spm_vec(par_rep{i})); r2loop(i)]';
            parBank = [parBank parI(:,i) ];
        end
    end
    [Ylist Ilist] = sort(r2loop,'descend');
    bestr2(ii) = Ylist(1);
    %     Ilist(isnan(r2loop(Ilist))) = []; % reconstruct Ilist without NaNs
    
    % Clip parBank to the best (keeps size manageable
    [dum V] = sort(parBank(end,:),'descend');
    if size(parBank,2)>2^13
        parBank = parBank(:,V(1:2^12));
    else
        parBank = parBank(:,V);
    end
    
    %% Find error threshold for temperature (epsilon)
    parOptBank = parBank(:,parBank(end,:)>eps_exp);
    A = parOptBank(pIndMap,:);
    B = eig(cov(A));
    C = B/sum(B);
    fprintf('effective rank of optbank is %.0f\n',sum(cumsum(C)>0.01))
    if size(parOptBank,2)> R.SimAn.minRank-1
        if size(parOptBank,2) < 4*(R.SimAn.minRank-1)
            disp('Bank satisfies current eps')
            eps_act = eps_exp;
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        else % if the bank is very large than take subset
            disp('Bank is large taking new subset to form eps')
            parOptBank = parBank(:,intersect(1:4*R.SimAn.minRank,1:size(parBank,2)));
            eps_act = min(parOptBank(end,:));
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        end
    elseif itry < 3
        fprintf('Trying for the %.0f\n time with the current eps \n',itry)
        disp('Trying once more with current eps')
        if isfield(Mfit,'Rho')
            cflag = 1;
        end
        itry = itry + 1;
    elseif itry >= 3
        disp('Recomputing eps from parbank')
        parOptBank = parBank(:,intersect(1:2*R.SimAn.minRank,1:size(parBank,2)));
        eps_act = min(parOptBank(end,:));
        cflag = 1;
        itry = 0;
    end
    if itry==0
        % Compute expected gradient for next run
        delta_exp = eps_exp-eps_prior;
        fprintf('Expected gradient was %0.2f \n',delta_exp)
        delta_act = eps_act-eps_prior;
        fprintf('Actual gradient was %0.2f \n',delta_exp)
        eps_exp = eps_act + delta_act;
        fprintf('Exp-Act gradient was %0.2f \n',delta_exp-delta_act)
        % Save eps history and make actual eps new prior eps
        eps_prior = eps_act;
    end
    eps_rec(ii) = eps_act;
    %% Compute Posterior
    if cflag == 1 && itry == 0 % estimate new copula
        [Mfit,cflag] = postEstCopula(R,parOptBank,pInd,pIndMap,pOrg);
        par = postDrawCopula(Mfit,pOrg,pIndMap,pSigMap,rep);
    elseif cflag == 1 && itry <= 3 % Draw from old copula
        par = postDrawCopula(Mfit,pOrg,pIndMap,pSigMap,rep);
    else % Draw from Normal Distribution
        if ~isempty(parOptBank)
            Mfit.Mu = mean(parOptBank(pMuMap,:),2);
            Mfit.Sigma = cov(parOptBank(pMuMap,:)');
        else
            Mfit.Mu = mean(parBank(pMuMap,intersect(1:2*R.SimAn.minRank,1:size(parBank,2))),2);
            Mfit.Sigma = cov(parBank(pMuMap,intersect(1:2*R.SimAn.minRank,1:size(parBank,2)))');
        end
        par = postDrawMVN(Mfit,pOrg,pIndMap,pSigMap,rep);
    end
    kldHist(ii)= KLDiv(R,p,m,parOptBank);

    parHist(ii) = averageCell(par);
    saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'],parHist)
    banksave{ii} = parBank(end,parBank(end,:)>eps_act);
    saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'],banksave)
    %%%%%%%%%%%%%%% SAVE PROGRESS, PLOTTING ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Ilist,2)>2 && R.plot.flag ==1
        if isfield(R.plot,'outFeatFx')
            %% Plot Data Features Outputs
            figure(1)
            clf
            fx = R.plot.outFeatFx;
            if size(Ilist,2)<12; xn = size(Ilist,2); else; xn = 12; end
            try
                fx({R.data.feat_emp},{feat_sim_rep{Ilist(1:xn)}},R.data.feat_xscale,R,1,[])
                drawnow; shg
            end
        end
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Plot parameters changes and tracking of fit
        if isfield(Mfit,'Rho')
            pmean = Mfit.Pfit;
        else
            pmean = p;
        end
        
        figure(2);    clf
        optProgPlot(1:ii,r2loop(Ilist(1)),pmean,banksave,eps_rec,bestr2,pInd,pSig,R)
        drawnow;shg
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Plot example time series
        figure(4)
        tvec_obs = R.IntP.tvec;
        tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
        R.IntP.tvec_obs = tvec_obs;
        ptr(1) = subplot(2,1,1);
        try
            plot(repmat(R.IntP.tvec_obs,size(xsims_rep{Ilist(1)}{1},1),1)',xsims_rep{Ilist(1)}{1}');
            xlabel('Time (s)'); ylabel('Amplitude')
            if numel(xsims_rep{Ilist(1)})>1
                ptr(2) = subplot(2,1,2);
                plot(repmat(R.IntP.tvec_obs,size(xsims_rep{Ilist(1)}{2},1),1)',xsims_rep{Ilist(1)}{2}'); %xlim([15 20])
                linkaxes(ptr,'x'); %xlim([10 20])
            else
                ptr(2) = subplot(2,1,2);
                plot(repmat(R.IntP.tvec_obs,size(xsims_rep{Ilist(1)}{1},1),1)',xsims_rep{Ilist(1)}{1}');
                xlim  ([8 10])
            end
            xlabel('Time (s)'); ylabel('Amplitude')
            legend(R.chsim_name)
            drawnow;shg
        end
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Export Plots
        if isequal(R.plot.save,'True')
            saveSimAnFigures(R,ii)
        end
    end
    disp({['Current R2: ' num2str(r2loop(Ilist(1)))];[' Temperature ' num2str(ii) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
    %% Save data
    if rem(ii,10) == 0
        saveSimABCOutputs(R,Mfit,m,parBank)
    end
    % Or to workspace
    %     assignin('base','R_out',R)
    
    if delta_act < R.SimAn.convIt
        disp('Itry Exceeded: Convergence')
        saveSimABCOutputs(R,Mfit,m,parBank)
        if R.plot.flag == 1
            saveSimAnFigures(R,ii)
        end
        return
    end
    
    
    ii = ii + 1;
    %     uv = whos;
    %     saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\memDebug_' R.out.tag '_' R.out.dag '.mat'],uv)
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
end
