function [R,parBank] = SimAn_ABC_211218(ic,~,p,m,R,parBank)
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
% workspace. The field R.mfit is appended with the Rho (normalized covariance) and nu
% (degrees of freedom) of the estimated copula.
% There are several plotting functions which will track the progress of the
% annealing.
%
% Timothy West (2018) - UCL CoMPLEX
% / UCL, Wellcome Trust Centre for Human Neuroscience
%%%%%%%%%%%%%%%%%%%%%%
%% Setup for annealing
if nargin<6
    parBank = [];
end
if isempty(m)
    m.m = 1;
end
pOrg = p; % Record prior parameters.

% Set Fixed Annealing Parameters
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
ii = 1;
eps_prior = -3; % prior eps (needed for gradient approximation);
eps_exp = -2; eps_act = eps_prior;
cflag = 0; delta_act = 1; 
% Compute indices of parameters to be optimized
[pInd,parMu,parSig] = parOptInds_110817(R,p,m.m); % in structure form
% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*3); %Ensure rank of sample is large enough to compute copula
% set initial batch of parameters from gaussian priors
if isfield(R,'Mfit')
    rep = repset;
    par = postDrawCopula(R.Mfit,p,pIndMap,rep);
else
    rep = repset;
        Mfit.Mu = parMu(pIndMap);
        Mfit.Sigma = diag(parSig(pIndMap));
        par = postDrawMVN(Mfit,pOrg,pIndMap,rep);
end
itry = 0;
%% Main Annealing Loop
while ii <= searchN
    rep = repset;
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
        if size(parOptBank,2) < 2*(R.SimAn.minRank-1)
            disp('Bank satisfies current eps')
            eps_act = eps_exp;
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        else % if the bank is very large than take subset
            disp('Bank is large taking new subset to form eps')
            parOptBank = parBank(:,1:2*R.SimAn.minRank);
            eps_act = mean(parOptBank(end,:));
            cflag = 1; % copula flag (enough samples)
            itry = 0;  % set counter to 0
        end
    elseif itry < 2
        disp('Trying once more with current eps')
        if isfield(Mfit,'Rho')
            cflag = 1;
        end
        itry = itry + 1;
    elseif itry > 1
        disp('Recomputing eps from parbank')
        parOptBank = parBank(:,2:R.SimAn.minRank);
        eps_act = mean(parOptBank(end,:));
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
        par = postDrawCopula(Mfit,pOrg,pIndMap,rep);
    elseif cflag == 1 && itry<= 1 % Draw from old copula
        par = postDrawCopula(Mfit,pOrg,pIndMap,rep);
    else % Draw from Normal Distribution
        Mfit.Mu = mean(parBank(pIndMap,:),2);
        Mfit.Sigma = cov(parBank(pIndMap,:)');
        par = postDrawMVN(Mfit,pOrg,pIndMap,rep);
    end
    
    %%%%%%%%%%%%%%% SAVE PROGRESS, PLOTTING ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Ilist,2)>2
        np = round(0.15*size(Ilist,2));
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
        banksave{ii} = parBank(:,parBank(end,:)>eps_act);
        figure(2);    clf
        optProgPlot(1:ii,r2loop(Ilist(1)),pmean,banksave,eps_rec,bestr2,pInd,R)
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
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\feattrack\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\r2track\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\dist_track\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
            %     close all
            pathstr = [R.rootn 'outputs\' R.out.tag '\' R.out.dag '\OptSaves\'];
            if ~exist(pathstr, 'dir')
                mkdir(pathstr);
            end
            disp({['Current R2: ' num2str(r2loop(Ilist(1)))];[' Temperature ' num2str(ii) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
        end
    end
    %     disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
    %% Save data
    if rem(ii,10) == 0
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'],Mfit)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'],m)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'],R)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'],parBank)
    end
    % Or to workspace
    %     assignin('base','R_out',R)
    
    if delta_act < 5e-3
        disp('Itry Exceeded: Convergence')
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'],Mfit)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'],m)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'],R)
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'],parBank)
        return
    end
    
    
    ii = ii + 1;
    %     uv = whos;
    %     saveMkPath([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\memDebug_' R.out.tag '_' R.out.dag '.mat'],uv)
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
end
