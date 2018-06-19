clear all; close all
rng(2017)
load bgc_params
load('C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\rat_121016\data\processed\L22_lesion_rat_121016.mat')

chloc_name = {'STR','GP','STN'};
data = (1).*FTdata.ContData.trial{1}; fsamp = FTdata.fsample;

[F_data,meancsd_data] = condenseConMat(data,chloc_name,FTdata.label',fsamp,9,6:0.2:30);
save F_data
 load F_data
M.dipfit.model(1) = [];
M.dipfit.model(2:10) = M.dipfit.model(1);
% [f csdtar] = synthspec_250117(18,18);
% csdtar = csdtar.*1.2e-5;
% csdtar = normaliseV(csdtar);
% fx_bgc(x,u,P,M)
% x = states
% u = input
% P = parameters
% M = model (fixed parameters)

%% Simulation length
tend = 15;
fsamp = 1000;
brn = .5; % burn in (seconds)
N = tend*fsamp;
tvec = linspace(0,tend,N);
rtol = 0.9;
% States of interest (membrane voltage)
v_inds = [1 3 5 7 9];

%% Inputs
U.dt = 1/fsamp;
U.u  = randn(N,5)/8;   % WHITE
% U.u = (3+colornoise(0,1,tend,N)/2).*U.u;
% U.u  = sqrt(spm_Q(1/16,N))*U.u; % PINK
% U.u  = (1/128)*sin(2*pi*15*tvec); % SINE
% U.u = U.u';



% Uex.u  = (sqrt(spm_Q(1/16,N))*Uex.u); % PINK
% Uex.u  = (1/128)*sin(2*pi*18*tvec);
% Uex.u = Uex.u';
M.Uen_C = ones(1,5);
M.LFPgain = [2 2 2 2 2];

P.Uen_C = zeros(1,5);
P.T = zeros(1,5);
P.G = zeros(1,5);
P.LFPgain = repmat(0,1,5);
%% START MAIN LOOP
searchN = 100;
Tm = 1; % temperature
alpha = 0.99; % alpha increment
rep = 24;
for ii = 1:searchN
    clear r2loop
    stdev = (1/1.5)*Tm(ii); % Set maximum width of distribution from which to mutate
    for jj = 1:rep % Replicates for each temperature
        Pnew = P;
                % Melt Parameters
        if ii>0
            jc = jj;
            xMax = 4; xMin = -4; % Set maximum distance moved from prior
            
            T = sparse(normrnd(P.T,stdev,1,size(P.T,2)));
            T(T>xMax) = xMax; T(T<xMin) = xMin;
            Pnew.T = T;
            
            G = sparse(normrnd(P.G,stdev,1,size(P.G,2)));
            G(G>xMax) = xMax; G(G<xMin) = xMin;
            Pnew.G = G;
            
            Uen_C = normrnd(P.Uen_C,stdev/2,1,size(P.Uen_C,2));
            Uen_C(Uen_C>xMax) = xMax; Uen_C(Uen_C<xMin) = xMin;
            Pnew.Uen_C = Uen_C;
            
            LFPgain = normrnd(P.LFPgain,stdev,1,size(P.Uen_C,2));
            LFPgain(LFPgain>(xMax/2)) = xMax; LFPgain(LFPgain<xMin) = xMin;
            Pnew.LFPgain = LFPgain;
        end   
       
        % Initialise integrator
        tme = 0;
        x1 = full(M.x);
        xstore = zeros(size(v_inds,2),1);
        LFP = zeros(5,N-1);
        
        % Integrate (Heuns Method)
        for tstep = 1:N-1;
            fx = fx_bgc_hpd_enU(x1,U.u(tstep,:),Pnew,M);
%             fhan   = fcnchk('fx_bgc_hpd_enU','x','u','P','M');
%             %         fhan = inline(char('fx_bgc_hpd_enU')); %,'x','U','P','M');
%             D=1;
%             dfdx = spm_cat(spm_diff(fhan,x1,U.u(tstep,:),Pnew,M,1));
%             x1 = x1 + spm_dx(D*dfdx,D*fx,U.dt);
            x1 = x1 + (U.dt.*fx);
            LFP(:,tstep) = x1(:,v_inds);
            tme(tstep+1) = tme(tstep) + U.dt;
            disp(tme(tstep))
        end
        % Delete burn in
        LFP(:,1:round(brn*fsamp)) = [];
        tme(:,1:round(brn*fsamp)) = [];
        gain = M.LFPgain.*exp(Pnew.LFPgain);
        LFP = bsxfun(@times,LFP,gain');
        % Compute gain
        
        % compute power spectra
        [F_sim,meancsd_sim] = condenseConMat(LFP(1:3,:),chloc_name,{'STR','GPe','STN'},fsamp,10,F_data);
        r2loop(jj) =compareCSD(meancsd_data,meancsd_sim);
        
        pxxsav{jj} = meancsd_sim; fsav{jj} = F_sim;
        ffxsav{jj} = F_data;
        
        Pbest(jj) = Pnew;
    end
    
    % Bug fix - for bad simulations make sure not compared
    r2loop(isnan(r2loop)) = 0;
    r2loop(imag(r2loop)~=0) = 0;
    
    % Find best fitting replicate
    [Y I] = max(r2loop);
    tbr2(ii) = Y; % Save R2 value
    % Recompute temperature
    Tm(ii+1) = Tm(ii)*alpha;
    if ii>1
        % Compute temperature and acceptance probability
        pd = (tbr2(ii)-tbr2(ii-1))/tbr2(ii); % Rescaled difference between mutants
        %         a = exp(-pd/Tm(ii)); % acceptance probability
        a = exp(pd/Tm(ii)^1.2)/4;
        rej = rand>a;
        
        if tbr2(ii) == 0
            rej = 1;
        end
        if rej~=1
            P = Pbest(I);
        else
            %             tbr2(ii) = tbr2(ii-1)
        end
        if tbr2(ii)>0.9
            rep = 16;
        end
        if tbr2(ii)>rtol
            break
        end
    end
    close all
    % Plot results
        figure(1); clf
    connection_diagram(M,P);
    set(gcf,'Position',[626    99   506   611]); drawnow; shg
    saveallfiguresFIL(['C:\Users\twest\Documents\Work\PhD\LitvakProject\Bernadette_DCM\convert2ODE\searchgif\ratCSD\L22\bgc_simanneal_condiag_rat_' num2str(ii) '_'],'-jpg',1,'-r900'); close all

    figure(2); clf
    csdplotter(meancsd_data,pxxsav{I},F_data)
    %     plot(fsav(I,:)',pxxsav(I,:)','k','linewidth',3); hold on; plot(fsav(I,:)',ffxsav(I,:)','k--','linewidth',2); xlim([1 100]); ylim([0 1]);
    %     hold on;  plot(fsav(:,:)',pxxsav(:,:)')
    % %     if max(max(pxxsav(I,:)))>0.1; ylim([0 1.1*max(max(pxxsav(I,:)))]); end
    %     xlabel('frequency'); ylabel('power');
    title(['Temperature: ' num2str(Tm(ii)) ' R2:' num2str(tbr2(ii))])
    set(gcf,'Position',[824   137   918   690])
    drawnow; shg
    saveallfiguresFIL(['C:\Users\twest\Documents\Work\PhD\LitvakProject\Bernadette_DCM\convert2ODE\searchgif\ratCSD\L22\bgc_simanneal_powspec_rat_' num2str(ii) '_'],'-jpg',1,'-r900'); close all
    
    clear fsav pxxsav
end
% Plot anneal analysis
figure(3)
plot(Tm(1:end-1),tbr2); shg
Pfit = P;
save('C:\Users\twest\Documents\Work\PhD\LitvakProject\Bernadette_DCM\convert2ODE\searchgif\ratCSD\L22\Pfitted_simanneal_ratfit','Pfit','M','U')