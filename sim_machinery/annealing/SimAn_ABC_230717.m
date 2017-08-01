function [R] = SimAn_ABC_230717(ic,u,p,m,R)
% This annealing function uses approximate Bayesian computation in order to
% estimate the posterior parameter distributions, using a shifting epsilon
% that moves with temperature
if isempty(m)
    m.m = 1;
end
parBank = []; parOptBank = []; eps = -1;  iflag = 0; par = cell(1,R.SimAn.rep); psave = []; itry = 0;
%% Notes
% This version collapses parameter resamples into a single function
% adapted for generic usage
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;
alpha = R.SimAn.alpha;
ii = 1;
while ii <= searchN
    stdev = (R.SimAn.jitter*Tm(ii)); % Set maximum width of distribution from which to mutate
    rep = repset;
    parfor jj = 1:rep % Replicates for each temperature
        x = ic;
        %% Resample Parameters
        if iflag == 1
            pnew = par{jj};
        else %ifjj > 1
            pnew = resampleParameters_240717(R,p,stdev,m.m); % Draw from prior
            %         else
            %             pnew = p;
        end % only do outside first run
         
        %% Simulate New Data
        % Integrate in time master fx function
        %         xsims = eval([R.IntP.intFx R.IntP.intFxArgs]);
        xsims = R.IntP.intFx(x,m,pnew,R,u);
        
        if isfield(R.obs,'obsFx') % Run Observer function
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        if isfield(R.obs,'transFx') % Run Data Transform
            %% Construct CSD and compare to data
            %             fx = R.obs.transFx;
            %             [~,feat_sim] = fx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,10,R);
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        % Now using NRMSE
        r2mean  = R.IntP.compFx(R,feat_sim);
%         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        feat_sim_rep{jj} = feat_sim;
        disp(['Iterate ' num2str(jj) ' temperature ' num2str(ii)])
    end % Iterates per temperature
    
    
    %% Evaluate Results and Compute Acceptance Probability
    % Retrieve fits
    r2loop = [r2rep{:}];
    %5 Make failed simulations
    r2loop(r2loop==1) = -1e3;
    r2loop(isnan(r2loop)) = -1e3;
    r2loop(imag(r2loop)~=0) = -1e3;
    % Find best fitting replicate
    %     [Y I] = max(r2loop);
    [Ylist Ilist] = sort(r2loop,'descend');
    I = Ilist(1);
    clear pmean
    p = par_rep{I};
    psave{ii} = p;
    tbr2(ii) = Ylist(1);
    
    % Create bank of params and fits
    for i = 1:numel(par_rep)
        parBank = [parBank [full(spm_vec(par_rep{i})); r2loop(i)]];
    end
    
    % Find error size (epsilon)
    if itry<5
        if size(parBank,2)>R.SimAn.minRank
            eps = prctile(parBank(end,end-(rep/2):end),70);
            eps_tmp = (-3.*Tm(ii))+2;
            if eps>-0.1 && eps<eps_tmp
                eps = eps_tmp;
            end
            %+(Tm(ii)/5);
            eps = max(eps);
            parOptBank = parBank(:,parBank(end,:)>eps);
            while size(parOptBank,2)<R.SimAn.minRank
                eps = eps-0.5;
                parOptBank = parBank(:,parBank(end,:)>eps);
                if eps<-100
                    parOptBank = [];
                    break
                end
            end
        end
        assignin('base','parOptBank',parOptBank)
        assignin('base','parBank',parBank)
        disp(size(parOptBank,2))
        
        if ~isempty(parOptBank)
            j = 0;
            clear copU xf ilist
            for i = 1:(size(parOptBank,1)-1)
                x = parOptBank(i,:);
                if sum(abs(diff(x)))>0 && var(x)>1/8
                    j = j +1;
                    ilist(j) = i;
                    copU(j,:) = ksdensity(x,x,'function','cdf');
                    xf(j,:) = x;
                end
            end
        end
        if size(parOptBank,2)> R.SimAn.minRank-1 % Ensure rank is larger than parameters
            itry = 0;
            %         if size(parOptBank,2)>16
            %             parOptBank = parOptBank(end-24:end,:);
            %         end
            iflag = 1;
            [Rho,nu] = copulafit('t',copU','Method','ApproximateML');
            R.Mfit.Rho = Rho;
            R.Mfit.nu = nu;
            R.Mfit.tbr2 = tbr2;
            
            figure(3)
            clf
            plotDistChange_KS(Rho,nu,xf,psave{1},R)
            % If redraw then increase epsilon
            Tm(ii+1) = Tm(ii)*alpha;
        else % If not enough samples - redraw from old distribution
            Tm(ii+1) = Tm(ii);
            %          searchN = searchN;
            itry = itry +1;
            disp(['Trying to sample for epsilon ' num2str(eps) ', but cannot generate enough samples. Redrawing from distribution for ' num2str(itry) ' time'])
        end
    elseif itry < 10
        itry = 0;
        Tm(ii+1) = Tm(ii-2);
    else
        disp(['Stuck in minima, aborting :`('])
        break
    end
    
    if ii>1 || iflag == 1
        r = copularnd('t',Rho,nu,rep);
        clear x1
        for Q = 1:size(xf,1)
            x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
        end
        clear base
        base = spm_vec(p);
        for i = 1:rep
            base(ilist) = x1(:,i);
            par{i} = spm_unvec(base,p);
        end
    end
    
    
    
    %% Plot Feature Outputs
    figure(1)
    clf
    fx = R.plot.outFeatFx;
    if ii<2
        fx({R.data.feat_emp},{feat_sim_rep{I}},R.data.feat_xscale,R,1)
    else
        fx({R.data.feat_emp},{feat_sim_rep{Ilist(1:4)}},R.data.feat_xscale,R,I)
    end
    drawnow; shg
    
    % Plot parameters and R2 track
    figure(2)
    clf
    optProgPlot(Tm(1:ii),tbr2(1:ii),psave{ii},R)
    drawnow;shg
    %
    if istrue(R.plot.save)
        saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\feattrack\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
        saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\r2track\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
        saveallfiguresFIL_n([R.rootn 'outputs\csd_gif\dist_track\' sprintf('%d',[R.d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
        %     close all
        saveMkPath([R.rootn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],R)
        saveMkPath([R.rootn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],parBank)
        
    end
    disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag; ['Eps ' num2str(eps)]})
    
    %     if rem(ii,10) == 0
    %             saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],R)
    %             saveMkPath([R.rootn '\' R.projectn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' sprintf('%d',[R.d(1:3)]) '.mat'],parBank)
    %     end
    
    if iflag == 1
        ii = ii + 1;
    end
end

%     end
%% Convergence criteria
%     if ii>5
%         grad = max(tbr2(ii-5:ii)) - min(tbr2(ii-5:ii));
%         disp(sprintf('Optimization gradient over 5 steps is %d',grad))
%         if grad< R.SimAn.gradtol(1)
%             repset = 8;
%         elseif grad < R.SimAn.gradtol(2)
%             disp('Model Reached Convergence')
%             return
%         end
%     end