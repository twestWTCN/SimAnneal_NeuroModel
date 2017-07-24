function [R] = SimAn_ABC_230717(ic,u,p,m,R)
% This annealing function uses approximate Bayesian computation in order to
% estimate the posterior parameter distributions, using a shifting epsilon
% that moves with temperature
if isempty(m)
    m.m = 1;
end
parBank = []; parOptBank = []; eps = -1;  iflag = 0; par = cell(1,32);
%% Notes
% This version collapses parameter resamples into a single function
% adapted for generic usage
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;

pSkewSave{1} = spm_unvec(zeros(size(spm_vec(p))),p);
pPrecSave{1} = spm_unvec(ones(size(spm_vec(p))),p);
for ii = 1:searchN
    stdev = (R.SimAn.jitter*Tm(ii)); % Set maximum width of distribution from which to mutate
    if ii == 1
        rep = 1;
    else
        rep = repset;
    end
    
    if ii>1 % Avoids broadcasting whole history to the cluster
        PSS = pSkewSave{ii-1};
        PPS = pPrecSave{ii-1};
    else
        PSS = pSkewSave{ii};
        PPS = pPrecSave{ii};
    end
    parfor jj = 1:rep % Replicates for each temperature
        x = ic;
        %% Resample Parameters
        if iflag == 1
            pnew = par{jj};
        elseif jj > 1
            pnew = resampleParameters_050717(R,p,PSS,PPS,stdev,m.m); % Resample parameters
        else
            pnew = p;
        end % only do outside first run
        
        %% Simulate New Data
        % Integrate in time master fx function
        %         xsims = eval([R.IntP.intFx R.IntP.intFxArgs]);
        xsims = R.IntP.intFx(x,m,pnew,R);
        
        if isfield(R.obs,'obsFx') % Run Observer function
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        if isfield(R.obs,'transFx') % Run Data Transform
            %% Construct CSD and compare to data
            %             fx = R.obs.transFx;
            %             [~,feat_sim] = fx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,10,R);
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,7,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        % Now using NRMSE
        r2mean  = R.IntP.compFx(R,feat_sim);
        R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        feat_sim_rep{jj} = feat_sim;
        disp(['Iterate ' num2str(jj) ' temperature ' num2str(ii)])
    end % Iterates per temperature
    
    
    %% Evaluate Results and Compute Acceptance Probability
    % Retrieve fits
    r2loop = [r2rep{:}];
    %5 Make failed simulations
    r2loop(isnan(r2loop)) = -32;
    r2loop(imag(r2loop)~=0) = -32;
    
    % Create bank of params and fits
    if ii>2
        for i = 1:numel(par_rep)
            parBank = [parBank [spm_vec(par_rep{i}); r2loop(i)]];
            eps = prctile(parBank(end,:),60)+(Tm(ii)/20);
            parOptBank = parBank(:,parBank(end,:)>eps);
            
        end
    end
    
    if size(parOptBank,2)>16
        if size(parOptBank,2)>24
            parOptBank = parOptBank(end-24:end,:);
        end
        iflag = 1;
        j = 0;
        clear u xf
        for i = 1:(size(parOptBank,1)-1)
            x = parOptBank(i,:);
            
            if sum(abs(diff(x)))>0
                j = j +1;
                ilist(j) = i;
                u(j,:) = ksdensity(x,x,'function','cdf');
                xf(j,:) = x;
            end
        end
        [Rho,nu] = copulafit('t',u','Method','ApproximateML');
        figure(3)
        clf
        plotDistChange_KS(Rho,nu,xf,par_rep{1},pPrecSave,pSkewSave,stdev,R)
        
        r = copularnd('t',Rho,nu,rep);
        clear x1
        for Q = 1:size(xf,1)
            x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
        end
        clear base
        base = spm_vec(p);
        for i = 1:rep
            base(ilist) = x1(:,i)
            par{i} = spm_unvec(base,p);
        end
    end
    % Find best fitting replicate
    %     [Y I] = max(r2loop);
    [Ylist Ilist] = sort(r2loop,'descend');
    I = Ilist(1);
    alpha = R.SimAn.alpha;
    clear pmean
    
    psave(ii) = par_rep{I};
    tbr2(ii) = Ylist(1);
    %     end
    
    if ii>1
        % Compute temperature and acceptance probability
        %         pd = tbr2(ii)-tbr2(ii-1); %
        pd = (tbr2(ii)-tbr2(ii-1))/tbr2(ii); % Rescaled difference between mutants
        [pPrecSave,pSkewSave] = skewResamp(tbr2,psave,pPrecSave,pSkewSave,ii,R);
        %                  a = exp(-pd/Tm(ii)); % acceptance probability
        a = exp(pd/Tm(ii)^1.2)/2;
        acpt = a>rand;
        
        if tbr2(ii) == 0
            acpt = 0;
        end
        if acpt==1
            p = par_rep{I};
            disp('Parameters Accepted')
        else
            [y i] = max(tbr2);
            p = psave(i);
        end
        %         if tbr2(ii)>R.SimAn.rtol_repeat
        %             repset = R.SimAn.rep(2);
        %         end
        %         R = fitDepSimPar(R,tbr2(ii)); % Set fit dependent sim parameters
        if (tbr2(ii)/m.m)>R.SimAn.rtol_converge
            disp('Model Reached Convergence')
            break
        end
        %         figure(3)
        %         clf
        %         R.plot.distchangeFunc(R,psave,pPrecSave,pSkewSave,stdev,ii)
    end
    
    if ii>5
        grad = max(tbr2(ii-5:ii)) - min(tbr2(ii-5:ii));
        disp(sprintf('Optimization gradient over 5 steps is %d',grad))
        if grad< R.SimAn.gradtol(1)
            repset = 8;
        elseif grad < R.SimAn.gradtol(2)
            disp('Model Reached Convergence')
            return
        end
    end
    Tm(ii+1) = Tm(ii)*alpha;
    d = clock;
    figure(1)
    clf
    fx = R.plot.outFeatFx;
    if ii<2
        fx({R.data.feat_emp},{feat_sim_rep{I}},R.data.feat_xscale,R,1)
    else
        fx({R.data.feat_emp},{feat_sim_rep{Ilist(1:4)}},R.data.feat_xscale,R,I)
    end
    drawnow; shg
    
    figure(2)
    clf
    optProgPlot(Tm(1:ii),tbr2,psave(ii),R)
    
    %     ylim([0 1])
    
    drawnow;shg
    %
    if istrue(R.plot.save)
        saveallfiguresFIL_n([R.rootn '\' R.projectn '\outputs\csd_gif\feattrack\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
        saveallfiguresFIL_n([R.rootn '\' R.projectn '\outputs\csd_gif\r2track\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
        saveallfiguresFIL_n([R.rootn '\' R.projectn '\outputs\csd_gif\dist_track\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
        %     close all
    end
    disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    
end

% R.out.CSD = csd_rep{I};
% R.out.P = p;

% Plot CSD
% csdplotter(R.data.meancsd_data,meancsd_sim,F_sim,R)
