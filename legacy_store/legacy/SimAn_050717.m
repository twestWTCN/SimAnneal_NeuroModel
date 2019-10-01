function [R] = SimAn_050717(ic,u,p,m,R)
%% Notes
% This version collapses parameter resamples into a single function

searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;

pSkewSave{1} = spm_unvec(zeros(size(spm_vec(p))),p);
pPrecSave{1} = spm_unvec(ones(size(spm_vec(p))),p);
for ii = 1:searchN
    stdev = (1)*((0.95*Tm(ii)^2)+0.01); % Set maximum width of distribution from which to mutate
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
        if ii>1
                pnew = resampleParameters_050717(R,p,PSS,PPS,stdev,m.m); % Resample parameters
        else
                pnew = p;
        end % only do outside first run
        
        %% Simulate New Data
        % Integrate in time master fx function
        xstore = stepintegrator_delay(R,x,u,m,pnew);
        % Run Observer function
        xsims = observe_data(xstore,m,pnew,R);
        %% Construct CSD and compare to data
        [F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,[],R);
        % Now using NRMSE
        r2mean = compareCSD_240517(R.data.meancsd_data,meancsd_sim,R);
        % Now using NRMSE
        %         csdplotter_220517([],{meancsd_sim},R.frqz,R)
        
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        csd_rep{jj} = meancsd_sim;
        disp(['Iterate ' num2str(jj) ' temperature ' num2str(ii)])
    end % Iterates per temperature
    
    
    %% Evaluate Results and Compute Acceptance Probability
    r2loop = [r2rep{:}];
    % Bug fix - for bad simulations make sure not compared
    r2loop(isnan(r2loop)) = 0;
    r2loop(imag(r2loop)~=0) = 0;
    
    % Find best fitting replicate
    %     [Y I] = max(r2loop);
    [Ylist Ilist] = sort(r2loop,'descend');
    I = Ilist(1);
    tbr2(ii) = Ylist(1); % Save R2 value
    % Recompute temperature
    alpha = R.SimAn.alpha;
    psave(ii) = par_rep{I};
    if ii>1
        % Compute temperature and acceptance probability
        %         pd = tbr2(ii)-tbr2(ii-1); %
        pd = (tbr2(ii)-tbr2(ii-1))/tbr2(ii); % Rescaled difference between mutants
        [pPrecSave,pSkewSave] = skewResamp(tbr2,psave,pPrecSave,pSkewSave,ii,R);
        %                  a = exp(-pd/Tm(ii)); % acceptance probability
        a = exp(pd/Tm(ii)^1.2)/4;
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
        figure(3)
        clf
        plotDistChange(R,psave,pPrecSave,pSkewSave,stdev,ii)
    end
    
    if ii>5
        grad = max(tbr2(ii-5:ii)) - min(tbr2(ii-5:ii));
        disp(sprintf('Optimization gradient over 5 steps is %d',grad))
        if grad< R.SimAn.gradtol(1)
            repset = 8;
        elseif grad < R.SimAn.gradtol(2)
            disp('Model Reached Convergence')
            break
        end
    end
    Tm(ii+1) = Tm(ii)*alpha;
    d = clock;
    figure(1)
    clf
    if ii<2
        csdplotter_220517({R.data.meancsd_data},{csd_rep{I}},R.frqz,R)
    else
        csdplotter_220517({R.data.meancsd_data},{csd_rep{Ilist(1:8)}},R.frqz,R,I)
    end
    drawnow; shg
    
    figure(2)
    clf
    subplot(2,1,1)
    plot(-Tm(1:ii),tbr2)
    subplot(2,1,2)
    par = full(spm_vec(p));
    par(par<-20) = 0;
    bar(par)
    ylim([1.2*min(par) 1.2*max(par)]);
    xlim([0 length(par)])
    xlabel('parameter')
    ylabel('Posterior')
    
    %     ylim([0 1])
    
    drawnow;shg
    % figure(3)
    saveallfiguresFIL_n(['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\csd_gif\csdtrack\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',1);
    saveallfiguresFIL_n(['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\csd_gif\r2track\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',2);
    saveallfiguresFIL_n(['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\csd_gif\dist_track\' sprintf('%d',[d(1:3)]) '\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_'],'-jpg',1,'-r100',3);
    %     close all
    disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']; R.out.tag})
    
end

R.out.CSD = csd_rep{I};
R.out.P = p;

% Plot CSD
% csdplotter(R.data.meancsd_data,meancsd_sim,F_sim,R)
