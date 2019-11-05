function [R] = SimAn_010617(ic,u,p,m,R)
searchN = R.SimAn.searchN;
repset =  R.SimAn.rep(1);
Tm(1) =  R.SimAn.Tm;
for ii = 1:searchN
    stdev = (1/R.SimAn.jitter)*Tm(ii); % Set maximum width of distribution from which to mutate
    if ii == 1
        rep = 1;
    else
        rep = repset;
    end
    for jj = 1:rep % Replicates for each temperature
        x = ic;
        %% Resample Parameters
        xMax = R.SimAn.maxdev; xMin = -R.SimAn.maxdev; % Set maximum distance moved from prior
        pnew = p;
        %         mnew = m;
        if ii>1
            for src = 1:length(m.m)
                
                pold = p;
                for opI = 1:numel(R.SimAn.opPar)
                    switch R.SimAn.opPar{opI}
                        case 'T' % Time constants
                            %             mold = m(src);
                            T_old = pold.int{src}.T;
                            T = normrnd(T_old,stdev,1,size(T_old,2));
                            T(T>xMax) = xMax; T(T<xMin) = xMin;
                            pnew.int{src}.T = T;
                            
                        case 'G' % Intrinsic synaptic gains
                            
                            G_old = pold.int{src}.G;
                            G = normrnd(G_old,stdev,1,size(G_old,2));
                            G(G>xMax) = xMax; G(G<xMin) = xMin;
                            pnew.int{src}.G = G;
                            
                        case 'C' % input gains
                            C = normrnd(pold.C,stdev,size(pold.C,1),1);
                            C(C>xMax) = xMax; C(C<xMin) = xMin;
                            pnew.C = C;
                            
                        case 'A' % extrinsic connections
                            
                            for Ai = 1:2
                                A_old = reshape(pold.A{Ai*2},1,[]);
                                A = normrnd(A_old,stdev,1,size(A_old,2));
                                A(A>(xMax)) = xMax; A(A<xMin) = xMin;
                                A = A.*reshape(m.Bmod{Ai},1,[]);
                                A(A==0) = -32;
                                A = reshape(A,m.m,m.m);
                                pnew.A{Ai*2} = A;
                                pnew.A{(Ai*2)-1} = A;
                            end
                            
                        case 'D' % Extrinsic delays (explicit)
                            
                            D_old = reshape(pold.D,1,[]);
                            D = normrnd(D_old,stdev,1,size(D_old,2));
                            D(D>(xMax)) = xMax; D(D<-2) = xMin; % changed so doesnt go too small
                            D = D.*reshape(m.Bmod{1}+m.Bmod{2},1,[]);
                            D(D==0) = -32;
                            D = reshape(D,m.m,m.m);
                            D(D>log(R.IntP.buffer/10)) = log(R.IntP.buffer/10); % Ensure not over buffer
                            pnew.D = D;
                            
                        case 'LF' % Lead field
                            LF_old = p.obs.LF;
                            LF = normrnd(LF_old,stdev/16,1,size(LF_old,2));
                            LF(LF>xMax) = xMax; LF(LF<xMin) = xMin;
                            pnew.obs.LF = LF;
                            
                        case 'mix' % Signal mixing
                            mix_old = pold.obs.mixing;
                            mix = normrnd(mix_old,stdev/4,1,size(mix_old,2));
                            mix(mix>xMax) = xMax; mix(mix<xMin) = xMin;
                            pnew.obs.mixing = mix;
                    end
                end
                %                 lfpgain = normrnd(pold.lfpgain,stdev,1,size(pold.lfpgain,2));
                %                 lfpgain(lfpgain>(xMax)) = xMax; lfpgain(lfpgain<xMin) = xMin;
                %                 pnew(src).lfpgain = lfpgain;
            end % resample parameters per source
        end % only do outside first run
        
        %% Simulate New Data
        % Integrate in time master fx function
        xstore = R.IntP.intFunc(R,x,u,m,pnew);
        % Run Observer function
        xsims = observe_data(xstore,m,pnew,R);
        %% Construct CSD and compare to data
        [F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.csd.pow2_sim,R);
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
    [Y I] = max(r2loop);
    tbr2(ii) = Y; % Save R2 value
    % Recompute temperature
    alpha = R.SimAn.alpha;
    psave(ii) = par_rep{I};
    if ii>1
        
        % Compute temperature and acceptance probability
        %         pd = tbr2(ii)-tbr2(ii-1); %
        pd = (tbr2(ii)-tbr2(ii-1))/tbr2(ii); % Rescaled difference between mutants
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
        if tbr2(ii)>R.SimAn.rtol_repeat
            repset = R.SimAn.rep(2);
        end
        
        %         if (tbr2(ii)/m.m)>R.SimAn.rtol_converge
        %             disp('Model Reached Convergence')
        %             break
        %         end
    end
    Tm(ii+1) = Tm(ii)*alpha;
    d = R.d;
    figure(1)
    clf
    csdplotter_220517({R.data.meancsd_data},{csd_rep{I}},R.frqz,R)
    %         csdplotter_220517({R.data.meancsd_data},csd_rep,R.frqz,R,I)
    drawnow; shg
    
    figure(2)
    clf
    subplot(2,1,1)
    plot(-Tm(1:ii),tbr2)
    figure(2)
    subplot(2,1,2)
    par = full(spm_vec(p));
    par(par<-20) = 0;
    bar(par)
    ylim([1.2*min(par) 1.2*max(par)]);
    xlim([0 length(par)])
    xlabel('parameter')
    ylabel('Posterior')
    
    if ii>R.SimAn.ntol
        grad = mean(diff(tbr2(ii-R.SimAn.ntol:ii)));
        if grad< R.SimAn.gradtol(1);
            disp('Optimisation Gradient is at Tolerance 1: reducing trials per temp')
        elseif grad< R.SimAn.gradtol(2);
            repset =  R.SimAn.rep(2);
            disp('Optimisation Gradient is at Tolerance 2: Finish!')
            break
        else
            disp(['Optimisation Gradient: ' sprintf('%d',grad)])
        end
    end
    %     ylim([0 1])
    drawnow;shg
    saveallfiguresFIL_n(['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\delayed\' sprintf('%d',[d(1:3)]) '\csdtrack\' R.out.tag '\bgc_siman_csdtrack_' num2str(ii) '_' R.out.tag],'-jpg',1,'-r100',1);
    saveallfiguresFIL_n(['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\delayed\' sprintf('%d',[d(1:3)]) '\r2track\' R.out.tag '\bgc_siman_r2track_' num2str(ii) '_' R.out.tag],'-jpg',1,'-r100',2);
    %     close all
    disp({['Current R2: ' num2str(tbr2(ii))];[' Temperature ' num2str(Tm(ii)) ' K']})
    R.out.CSD = csd_rep{I};
    R.out.P = p;
    assignin('base',R.SimAn.saveout,R)
end

% R.out.CSD = csd_rep{I};
% R.out.P = p;
% Plot CSD
% csdplotter(R.data.meancsd_data,meancsd_sim,F_sim,R)
