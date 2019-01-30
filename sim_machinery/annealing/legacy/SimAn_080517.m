function [xobs] = SimAn_080517(ic,u,p,m,R)
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;
for ii = 1:searchN
    stdev = (4)*Tm(ii); % Set maximum width of distribution from which to mutate
    if ii == 1
        rep = 1;
    else
        rep = repset;
    end
    for jj = 1:rep % Replicates for each temperature
        %% Resample Parameters
        xMax = 8; xMin = -8; % Set maximum distance moved from prior
        pnew = p;
        %         mnew = m;
        if ii>2
            for src = 1:length(m(1).nstates)
                pold = p(src);
                %             mold = m(src);
                T = normrnd(pold.T,stdev,1,size(pold.T,2));
                T(T>xMax) = xMax; T(T<xMin) = xMin;
                pnew(src).T = T;
                
                G = normrnd(pold.G,stdev,1,size(pold.G,2));
                G(G>xMax) = xMax; G(G<xMin) = xMin;
                pnew(src).G = G;
                
                C = normrnd(pold.C,stdev/2,1,size(pold.C,2));
                C(C>xMax) = xMax; C(C<xMin) = xMin;
                pnew(src).C = C;
                
                A = normrnd(pold.A,stdev,1,size(pold.A,2));
                A(A>(xMax)) = xMax; A(A<xMin) = xMin;
                pnew(src).A = A;
                
%                 lfpgain = normrnd(pold.lfpgain,stdev,1,size(pold.lfpgain,2));
%                 lfpgain(lfpgain>(xMax)) = xMax; lfpgain(lfpgain<xMin) = xMin;
%                 pnew(src).lfpgain = lfpgain;
            end % resample parameters per source
        end % only do outside first run
        
        %% Simulate New Data
        % Integrate in time master fx function
        xstore = zeros(6,R.IntP.nt);
        for t = 1:R.IntP.nt
            [x xobs] = simAn_master_fx_bgc(ic,u,pnew,m,R.IntP.dt,t);
            xstore(:,t) = xobs';
%             t
        end
        % Run Observer function
        xsims = observe_data(xstore,m,p,R);
        %% Construct CSD and compare to data
        [F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,9,R);
        r2mean = compareCSD_080517(R.data.meancsd_data,meancsd_sim);
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        csd_rep{jj} = meancsd_sim;
    end % Iterates per temperature
    
    
    %% Evaluate Results and Compute Acceptance Probability
    r2loop = [r2rep{:}];
    % Bug fix - for bad simulations make sure not compared
    r2loop(isnan(r2loop)) = 0;
    r2loop(imag(r2loop)~=0) = 0;
    
    % Find best fitting replicate
    [Y I] = max(r2loop);
    tbr2(ii) = Y; % Save R2 value
    % Recompute temperature
    alpha = R.SimAn.alpha;

    if ii>1
        % Compute temperature and acceptance probability
        pd = (tbr2(ii)-tbr2(ii-1))/tbr2(ii); % Rescaled difference between mutants
        %         a = exp(-pd/Tm(ii)); % acceptance probability
        a = exp(pd/Tm(ii)^1.2)/4;
        acpt = a>rand;
        
        if tbr2(ii) == 0
            acpt = 0;
        end
        if acpt==1
            p = par_rep{I};
             disp('Parameters Accepted')
        end
        if tbr2(ii)>R.SimAn.rtol_repeat
            repset = R.SimAn.rep(2);
        end
        if tbr2(ii)>R.SimAn.rtol_converge
            disp('Model Reached Convergence')
            break
        end
    end
    Tm(ii+1) = Tm(ii)*alpha;
    clf
    csdplotter(R.data.meancsd_data,csd_rep{I},R.frqz,R)
    drawnow; shg
    disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']})
end
% Plot CSD
% csdplotter(R.data.meancsd_data,meancsd_sim,F_sim,R)
