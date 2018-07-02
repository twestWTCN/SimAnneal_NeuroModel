function [R,parBank] = SimAn_ABC_110817(ic,u,p,m,R,parBank)
%%%% SIMULATED ANNEALING for APROXIMATE BAYESIAN COMPUTATION of NONLINEAR
%%%% DYNAMICAL MODELS
% ---- 11/08/17---------------------------
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
%
% This function will output R at each annealing loop (assigned into 'base'
% workspace. The field R.mfit is appended with the Rho (covariance) and nu
% (degrees of freedom) of the estimated copula.
% There are several plotting functions which will track the progress of the
% annealing.
%
%%%%%%%%%%%%%%%%%%%%%%
% Setup for annealing
if nargin<6
    parBank = [];
end
if isempty(m)
    m.m = 1;
end
pOrg = p; % Record prior parameters.

% Initialise variables
parOptBank = []; eps = -1;  iflag = 0; par = cell(1,R.SimAn.rep); psave = []; itry = 0;
% Set Fixed Annealing Parameters
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;
alpha = R.SimAn.alpha;
ii = 1;

% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*2);
% set initial batch of parameters from gaussian priors
stdev = R.SimAn.jitter*Tm(1); % set the global precision
if isfield(R,'Mfit')
    xf = R.Mfit.xf;
    rep = repset;
    disp('Drawing from copula...')
    r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,rep);
    clear x1
    for Q = 1:size(xf,1)
        x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
    end
    % setup pars from base
    clear base
    base = repmat(spm_vec(p),1,rep);
    for i = 1:rep
        base(pIndMap,i) = x1(:,i);
        par{i} = spm_unvec(base(:,i),p);
    end
    iflag = 1;
else
    for jj = 1:repset
        par{jj} = resampleParameters_240717(R,p,stdev,m.m); % Draw from prior
    end
end
itry = 0;
eps_p = -100;
%% Main Annealing Loop
while ii <= searchN
    rep = repset;
    %% Batch Loop for Replicates for Generation of Pseudodata
    % This is where the heavy work is done. This is run inside parfor. Any
    % optimization here is prime.
    
    parfor jj = 1:rep % Replicates for each temperature
        u = innovate_timeseries(R,m);
        u = u.*sqrt(R.IntP.dt);
        x = ic;
        %% Resample Parameters
        pnew = par{jj};
        %% Simulate New Data
        % Integrate in time master fx function
        try
        xsims = R.IntP.intFx(R,x,u,pnew,m);
        catch
            disp('Simulation failed')
            xsims = nan(1,3);
        end
        if sum(isnan(xsims(:))) == 0
            try
                % Run Observer function
                glorg = pnew.obs.LF;
                % Subloop is local optimization of the observer gain
                % Parfor requires parameter initialisation
                gainlist = 0; %linspace(-0.5,1,8);
                feat_sim = cell(1,length(gainlist));
                xsims_gl = cell(1,length(gainlist));
                r2mean = zeros(1,length(gainlist));
                for gl = 1:length(gainlist)
                    pnew.obs.LF = glorg + gainlist(gl);
                    if isfield(R.obs,'obsFx')
                        xsims_gl{gl} = R.obs.obsFx(xsims,m,pnew,R);
                    else
                        xsims_gl{gl} =xsims;
                    end
                    % Run Data Transform
                    if isfield(R.obs,'transFx')
                        [~,feat_sim{gl}] = R.obs.transFx(xsims_gl{gl},R.chloc_name,R.chsim_name,1/R.IntP.dt,R.obs.SimOrd,R);
                    else
                        feat_sim{gl} = xsims_gl{gl}; % else take raw time series
                    end
                    % Compare Pseudodata with Real
                    r2mean(gl)  = R.IntP.compFx(R,feat_sim{gl});
                end
                % F
                [r2 ir2] = max(r2mean);
                pnew.obs.LF = glorg + gainlist(ir2);
                %         disp(pnew.obs.LF)
                %         toc
                % plot if desired
                %                 R.plot.outFeatFx({R.data.feat_emp},{feat_sim{ir2}},R.data.feat_xscale,R,1,[])
            catch
                r2 = -inf;
                ir2 =1;
                xsims_gl{1} = NaN;
                feat_sim{1} = NaN;
            end
        else
            r2 = -inf;
            ir2 =1;
            xsims_gl{1} = NaN;
            feat_sim{1} = NaN;
        end
        
        r2rep{jj} = r2;
        par_rep{jj} = pnew;
        xsims_rep{jj} = xsims_gl{ir2};
        feat_sim_rep{jj} = feat_sim{ir2};
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
        if ~isnan(r2loop(i))
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
        parBank = parBank(:,V(1:2^13));
    else
        parBank = parBank(:,V);
    end
    
    % Find error threshold for temperature (epsilon)
    eps = NaN;
    if size(parBank,2)>R.SimAn.minRank
        %L = round(size(parBank,2)*0.4);
        %         eps_op(2) = prctile(parBank(end,:),85);
        %         eps = max(eps_op);
        %         try
        %             L = rep;
        %             eps = prctile(parBank(end,end-L:end),75); % percentile eps
        %         catch
        %             L = round(size(parBank,2)*0.4);
        %             eps = prctile(parBank(end,end-L:end),75); % percentile eps
        %         end
        eps_tmp(1) = (-3.*Tm(ii))+2; % temperature based epsilon (arbitrary function)
        clear parOptBank
        parOptBank = parBank(:,parBank(end,:)>eps_tmp(1));
        
        if size(parOptBank,2)< R.SimAn.minRank-1
            eps_check = [0 0];
            disp(['Annealing EPS not yielding large enough bank, try percentile EPS'])
            eps_tmp(2) = prctile(r2loop,85);
            A{2} = parBank(:,parBank(end,:)>eps_tmp(2));
            if size(A{2},2)>size(parOptBank,2)
                eps_check(1) = 1;
            end
            eps_tmp(3) = prctile(parBank(end,:),95);
            A{3} = parBank(:,parBank(end,:)>eps_tmp(3));
            if size(A{3},2)>size(parOptBank,2)
                eps_check(2) = 1;
            end
            Asel = [size(A{3},2)>size(A{2},2) size(A{3},2)<size(A{2},2)].*eps_check;
            dispmes = {'Using 85% ii Call Pars','Using 95% Parbank pars'};
            if any(Asel)
                parOptBank = A{find(Asel)+1};
                eps = eps_tmp(find(Asel)+1);
                disp(dispmes{find(Asel)})
            end
            disp(['90% EPS: ' num2str(eps) ', length ' num2str(size(parBank(:,parBank(end,:)>eps),2))])
        else
            itry = 0;
            eps = eps_tmp(1);
            disp(['Using Annealing EPS: ' num2str(eps) ', length ' num2str(size(parBank(:,parBank(end,:)>eps),2))])
        end
        clear A
        
        if  eps-eps_p < 0.005
            itry = itry + 1;
            if itry<=R.SimAn.copout(2)
                pip = 100;
                A = [1];
                while size(A,2)<R.SimAn.minRank-1
                    pip = pip-0.001;
                    eps = prctile(parBank(end,:),pip);
                    A = parBank(:,parBank(end,:)>eps);
                end
                parOptBank = A;
                clear A
            elseif (itry>R.SimAn.copout(2)) 
                pip = 100;
                A = [1];
                while size(A,2)<(R.SimAn.minRank)*0.60
                    pip = pip-0.001;
                    eps = prctile(parBank(end,:),pip);
                    A = parBank(:,parBank(end,:)>eps);
                end
                parOptBank = A;
                clear A
            end
        end
        
        eps_p = eps; % past value of eps (for next loop)
        
        % Form the bank of parameters exceeding acceptance level
        % If the size of the bank is less than the rank of the optimized
        % data and more than a given number of replicates are attempted
        % then choose epsilon to give bank exceeding rank (bare minimum for
        % copula formation
        if ((R.SimAn.minRank-size(parOptBank,2))< 0.35*R.SimAn.minRank) && itry>R.SimAn.copout(1)
            disp(['Rank is within limits, filling with pseudopars: ' num2str(eps)])
            aN = R.SimAn.minRank-size(parOptBank,2);
            mu = mean(parBank(:,parBank(end,:)>eps),2);
            sig = cov(parBank(:,parBank(end,:)>eps)');
            sig = (sig + sig.') / 2; % Ensures symmetry constraints
            try
                parOptBank = [parOptBank mvnrnd(mu,sig,aN)'];
            catch
                disp('covar matrix is non-symmetric postive definite')
            end
            itry =0;
        end
        
        % Crops the optbank to stop it getting to big (memory)
        if size(parOptBank,2)> 2^9
            [dum V] = sort(parOptBank(end,:),'descend');
            parOptBank = parOptBank(:,V(1:2^9));
        end
        
        % assign to base workspace in case stop early
        %         assignin('base','parOptBank',parOptBank)
        %         assignin('base','parBank',parBank)
        
        disp(['The number of parsets within epsilon is ' num2str(size(parOptBank,2))])
        
        %% This is where the multivariate posterior parameter estimate is computed using copulas
        if size(parOptBank,2)> R.SimAn.minRank-1 % Ensure rank is larger than parameters
            disp('Forming new copula...')
            clear copU xf ilist
            % First form kernel density estimates for each optimized
            % parameter
            clear copU
            for i = 1:size(pIndMap,1)
                x = parOptBank(pIndMap(i),:); % choose row of parameter values
                copU(i,:) = ksdensity(x,x,'function','cdf'); % KS density estimate per parameter
                xf(i,:) = x;
            end
            try
                [Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
                % Save outputs that specify the copula
                R.Mfit.Rho = Rho;
                R.Mfit.xf = xf;
                R.Mfit.nu = nu;
                R.Mfit.tbr2 = parOptBank(end,1); % get best fit
                R.Mfit.Pfit = spm_unvec(mean(parOptBank,2),p);
                Mfit_hist(ii) = R.Mfit;
                %%% Plot posterior, Rho, and example 2D/3D random draws from copulas
                figure(3)
                clf
                plotDistChange_KS(Rho,nu,xf,pOrg,pInd,R,stdev)
                %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
                 iflag = 1;
                % If redraw then increase epsilon
            catch
                disp('The estimate of Rho has become rank-deficient.  You may have too few data, or strong dependencies among variables.')
                    p = spm_unvec(mean(parOptBank,2),p);
                    iflag = 0;
            end
        else  % If not enough samples - redraw using means from the best fitting model
            %             p = par_rep{Ilist(1)}; % The best par_rep
            p = spm_unvec(mean(parOptBank,2),p);
            iflag = 0;
        end
    end
    eps_rec(ii) = eps;
    
    %% Now draw parameter sets for the next set of replicates
    if iflag == 0 && (itry < R.SimAn.copout(2))
        disp(['To many tries/trying old copula ' num2str(stdev)])
        try
            R.Mfit =  Mfit_hist(end);
            iflag = 1;
        catch
            disp(['No copula available, reverting to normal distribution on best fitting posteriors with stdev ' num2str(stdev)])
        end
    end
    if iflag == 1
        disp('Drawing from copula...')
        r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,rep);
        clear x1
        for Q = 1:size(xf,1)
            x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
        end
        % setup pars from base
        clear base
        base = repmat(spm_vec(p),1,rep);
        for i = 1:rep
            base(pIndMap,i) = x1(:,i);
            par{i} = spm_unvec(base(:,i),p);
        end
    elseif iflag == 0 && itry > R.SimAn.copout(2)
        disp(['No copula available, reverting to normal distribution on best fitting posteriors'])
        mu = mean(parOptBank(1:end-1,:),2);
        pm = spm_unvec(mu,p);
        sig = cov(parOptBank(1:end-1,:)');
        sig = (sig + sig.') / 2; % Ensures symmetry constraints
        stdev = mean(diag(sig));
        for jj = 1:repset
            par{jj} = resampleParameters_240717(R,pm,stdev,m.m); % Draw from prior
        end
        itry = 0;
    end
    
    %%%%%%%%%%%%%%% SAVE PROGRESS, PLOTTING ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Ilist,2)>2
        np = round(0.15*size(Ilist,2));
        %% Plot Data Features Outputs
        figure(1)
        clf
        fx = R.plot.outFeatFx;
        if size(Ilist,2)<12; xn = size(Ilist,2); else; xn = 12; end
        fx({R.data.feat_emp},{feat_sim_rep{Ilist(1:xn)}},R.data.feat_xscale,R,1,[])
        drawnow; shg
        
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Plot parameters changes and tracking of fit
        if exist('base')
            pmean = spm_unvec(mean(base,2),pOrg);
        else
            pmean = p;
        end
        banksave{ii} = parBank(:,parBank(end,:)>eps);
        figure(2);    clf
        optProgPlot(Tm(1:ii),r2loop(Ilist(1)),pmean,banksave,eps_rec,bestr2,pInd,R)
        drawnow;shg
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Plot example time series
        figure(4)
        tvec_obs = R.IntP.tvec;
        tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
        R.IntP.tvec_obs = tvec_obs;
        subplot(2,1,1)
        plot(repmat(R.IntP.tvec_obs,size(xsims_rep{Ilist(1)},1),1)',xsims_rep{Ilist(1)}');
        xlabel('Time (s)'); ylabel('Amplitude')
        subplot(2,1,2)
        plot(repmat(R.IntP.tvec_obs,size(xsims_rep{Ilist(1)},1),1)',xsims_rep{Ilist(1)}'); xlim([5 6])
        xlabel('Time (s)'); ylabel('Amplitude')
        legend(R.chsim_name)
        drawnow;shg
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Export Plots
        if istrue(R.plot.save)
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '\csd_gif\feattrack\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '\r2track\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
            saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '\dist_track\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
            %     close all
            pathstr = [R.rootn '\outputs\' R.out.tag '\' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '\OptSaves\'];
            if ~exist(pathstr, 'dir')
                mkdir(pathstr);
            end
            try
                save([pathstr 'modelfit_' R.out.tag '_' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '.mat'],'R','parBank','p','m','Mfit_hist')
            catch
                save([pathstr 'modelfit_' R.out.tag '_' sprintf('%d_%.2d_%.2d',[R.d(1:3)]) '.mat'],'R','parBank','p','m')
            end
            %             save([R.rootn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],'R')
            %             save([R.rootn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],'parBank')
            disp({['Current R2: ' num2str(r2loop(Ilist(1)))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
        end
    end
    %     disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
    %% Save data
        if rem(ii,10) == 0
                saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' R.out.dag '.mat'],R)
                saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' R.out.dag '.mat'],parBank)
        end
    % Or to workspace
    %     assignin('base','R_out',R)
    if iflag == 1
        Tm(ii+1) = Tm(ii)*alpha;
        ii = ii + 1;
    end
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
end
