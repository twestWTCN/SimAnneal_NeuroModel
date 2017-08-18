function [R] = SimAn_ABC_110817(ic,u,p,m,R)
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

if isempty(m)
    m.m = 1;
end

pOrg = p; % Record prior parameters.

% Initialise variables
parBank = []; parOptBank = []; eps = -1;  iflag = 0; par = cell(1,R.SimAn.rep); psave = []; itry = 0;

% Set Fixed Annealing Parameters
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;
alpha = R.SimAn.alpha;
ii = 1;

% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
% set initial batch of parameters from gaussian priors
stdev = R.SimAn.jitter*Tm(1); % set the global precision
for jj = 1:repset
    par{jj} = resampleParameters_240717(R,p,stdev,m.m); % Draw from prior
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
        u = u./R.IntP.dt;
        x = ic;
        %% Resample Parameters
        pnew = par{jj};
        %% Simulate New Data
        % Integrate in time master fx function
        xsims = R.IntP.intFx(R,x,u,pnew,m);
        % Run Observer function
        if isfield(R.obs,'obsFx')
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        % Run Data Transform
        if isfield(R.obs,'transFx')
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        % Compare Pseudodata with Real
        r2mean  = R.IntP.compFx(R,feat_sim);
        % plot if desired
        %         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        xsims_rep{jj} = xsims;
        feat_sim_rep{jj} = feat_sim;
        disp(['Iterate ' num2str(jj) ' temperature ' num2str(ii)])
    end % End of batch replicates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARAMETER OPTIMIZATION BEGINS HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Concatanate Batch Results and Decide Acceptance Level Epsilon
    % Retrieve fits
    r2loop = [r2rep{:}];
    % Delete failed simulations
    r2loop(r2loop==1) = NaN;
    r2loop(isnan(r2loop)) = NaN;
    r2loop(imag(r2loop)~=0) = NaN;
    % Append succesful replicates to bank of params and fits
    %(parameter table, with fits)
    for i = 1:numel(r2loop)
        if ~isnan(r2loop(i))
            parBank = [parBank [full(spm_vec(par_rep{i})); r2loop(i)]];
        end
    end
    [Ylist Ilist] = sort(r2loop,'descend');
    Ilist(isnan(r2loop(Ilist))) = []; % reconstruct Ilist without NaNs
    
    % Clip parBank to the best (keeps size manageable
    [dum V] = sort(parBank(end,:),'ascend');
    if size(parBank,2)>2048
        parBank = parBank(:,V(1:2048));
    else
        parBank = parBank(:,V);
    end
    
    % Find error threshold for temperature (epsilon)
    if size(parBank,2)>R.SimAn.minRank
        L = round(size(parBank,2)*0.4);
        eps = prctile(parBank(end,end-L:end),90); % percentile eps
        eps_tmp = (-3.*Tm(ii))+2; % temperature based epsilon (arbitrary function)
        if  eps-eps_p < 0.005
            eps = eps + 0.008;
            disp('EPS was getting stuck, forcing increase')
        end
        
        clear parOptBank
        parOptBank = parBank(:,parBank(end,:)>eps_tmp);
        if size(parOptBank,2)< R.SimAn.minRank-1
            disp(['Annealing EPS not yielding large enough bank, try percentile EPS'])
            A = parBank(:,parBank(end,:)>eps);
            if size(A,2)>size(parOptBank,2)
                parOptBank = A;
            end
        end
        eps_p = eps;
        disp(['90% EPS: ' num2str(eps) ', length ' num2str(size(parBank(:,parBank(end,:)>eps),2))])
        disp(['Anneal EPS: ' num2str(eps_tmp) ', length ' num2str(size(parBank(:,parBank(end,:)>eps_tmp),2))])
        % Form the bank of parameters exceeding acceptance level
        
        % If the size of the bank is less than the rank of the optimized
        % data and more than a given number of replicates are attempted
        % then choose epsilon to give bank exceeding rank (bare minimum for
        % copula formation
        if size(parOptBank,2)> R.SimAn.minRank-1
            itry = 0; % set itry to 0 as criterion is met
        elseif itry >=3
            disp(['To many attempts using old copula/dist, using pseudopars EPS instead: ' num2str(eps)])
            aN = R.SimAn.minRank-size(parOptBank,2);
            mu = mean(parBank(:,end-R.SimAn.minRank:end),2);
            sig = cov(parBank(:,end-R.SimAn.minRank:end)');
            sig = (sig + sig.') / 2; % Ensures symmetry constraints
            try
                parOptBank = [parOptBank mvnrnd(mu,sig,aN)'];
            catch
                disp('covar matrix is non-symmetric postive definite')
            end
            %             while size(parOptBank,2)<R.SimAn.minRank % ignore the epsilon if not enough rank
            %                 eps = eps-0.005;
            %                 parOptBank = parBank(:,parBank(end,:)>eps);
            %                 if eps<-100
            %                     parOptBank = [];
            %                     break
            %                 end
            %             end
            itry = 0;
        else
            itry = itry + 1;
            disp(['Difference between bank size and minimum rank is within tolerance, resampling for ' num2str(itry) ' time'])
        end
        
        % Crops the optbank to stop it getting to big (memory)
        if size(parOptBank,2)> 1024
            [dum V] = sort(parOptBank(end,:),'descend');
            parOptBank = parOptBank(:,V(1:1024));
        end
        
        % assign to base workspace in case stop early
        assignin('base','parOptBank',parOptBank)
        assignin('base','parBank',parBank)
        
        disp(['The number of parsets within epsilon is ' num2str(size(parOptBank,2))])
        
        %% This is where the multivariate posterior parameter estimate is computed using copulas
        if size(parOptBank,2)> R.SimAn.minRank-1 % Ensure rank is larger than parameters
            disp('Forming new copula...')
            iflag = 1; % flag to draw from distribution (should only go up once- following first copula)
            clear copU xf ilist
            % First form kernel density estimates for each optimized
            % parameter
            for i = 1:size(pIndMap,1)
                x = parOptBank(pIndMap(i),:); % choose row of parameter values
                copU(i,:) = ksdensity(x,x,'function','cdf'); % KS density estimate per parameter
                xf(i,:) = x;
            end
            [Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
            % Save outputs that specify the copula
            R.Mfit.Rho = Rho;
            R.Mfit.nu = nu;
            R.Mfit.tbr2 = parOptBank(end,1); % get best fit
            
            %%% Plot posterior, Rho, and example 2D/3D random draws from copulas
            figure(3)
            clf
            plotDistChange_KS(Rho,nu,xf,pOrg,pInd,R,stdev)
            %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
            
            % If redraw then increase epsilon
            iflag = 1;
        else  % If not enough samples - redraw using means from the best fitting model
            %             p = par_rep{Ilist(1)}; % The best par_rep
            p = spm_unvec(mean(parOptBank,2),p);
            
            disp(['Trying to sample for epsilon ' num2str(eps) ', but cannot generate enough samples. Redrawing from old copula for ' num2str(itry) ' time'])
        end
    end
    
    %% Now draw parameter sets for the next set of replicates
    if iflag == 1 && itry < 4
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
    else
        disp(['To many tries/ no copula available, reverting to normal distribution on best fitting posteriors with stdev ' num2str(stdev)])
        stdev = R.SimAn.jitter*Tm(ii); % set the global precision
        for jj = 1:repset
            par{jj} = resampleParameters_240717(R,p,stdev,m.m); % Draw from prior
        end
    end
    
    %%%%%%%%%%%%%%% SAVE PROGRESS, PLOTTING ETC. %%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Ilist,2)>2
        np = round(0.15*size(Ilist,2));
        %% Plot Data Features Outputs
        figure(1)
        clf
        fx = R.plot.outFeatFx;
        fx({R.data.feat_emp},{feat_sim_rep{Ilist(1:np)}},R.data.feat_xscale,R,1)
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
        optProgPlot(Tm(1:ii),r2loop(Ilist(1)),pmean,banksave,eps,pInd,R)
        drawnow;shg
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Plot example time series
        figure(4)
        tvec_obs = R.IntP.tvec;
        tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
        R.IntP.tvec_obs = tvec_obs;
        subplot(2,1,1)
        plot(repmat(R.IntP.tvec_obs,6,1)',xsims_rep{Ilist(1)}');
        xlabel('Time (s)'); ylabel('Amplitude')
        subplot(2,1,2)
        plot(repmat(R.IntP.tvec_obs,6,1)',xsims_rep{Ilist(1)}'); xlim([5 6])
        xlabel('Time (s)'); ylabel('Amplitude')
        legend(R.chsim_name)
        drawnow;shg
        %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
        %% Export Plots
        if istrue(R.plot.save)
            saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\feattrack\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
            saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\r2track\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
            saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\dist_track\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
            %     close all
            saveMkPath([R.rootn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],R)
            saveMkPath([R.rootn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],parBank)
            disp({['Current R2: ' num2str(r2loop(Ilist(1)))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
        end
    end
    %     disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
    %% Save data
    %     if rem(ii,10) == 0
    %             saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],R)
    %             saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],parBank)
    %     end
    % Or to workspace
    assignin('base','R_out',R)
   
Tm(ii+1) = Tm(ii)*alpha;
    if iflag == 1
        ii = ii + 1;
    end
    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%    %%%     %%%     %%%     %%%     %%%     %%%     %%%     %%%
end

