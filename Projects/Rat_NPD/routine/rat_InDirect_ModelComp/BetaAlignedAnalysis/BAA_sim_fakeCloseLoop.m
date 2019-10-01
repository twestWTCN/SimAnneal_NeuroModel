function [R] = BAA_sim_fakeCloseLoop(R,simtime,fresh)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,10,simtime);
P = permMod{1}.par_rep{1};


R.obs.csd.df = 0.25;
R = setSimTime(R,simtime);
cmap = brewermap(18,'Spectral');
PRC.condcmap = cmap; %([1 4 8 16 4 18],:);
PRC.condname = {'Fitted','1% STN->GPe','150% STN->GPe'}; %Fitted','1% M2->STN','150% M2->STN',

% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
fsamp = 1/R.IntP.dt;
conStren = [1]; %  0.001 1.3];
ck_1 = [];
% conStren =linspace(0.001,1.3,18);

if fresh
    
    phaseShift = linspace(0,2.*pi,18);
    %     intpow = []; maxpow = [];
    for p = 1:numel(phaseShift) %[1 10] %
        % Setup pulses for PRC computation
        uc_ip = {}; feat_sim = {}; xsim_ip = {};
        uc_ip{1} = uc;
        
        % Setup the simulations
        Pbase = P;
        %         Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*conStren(cond)); %
        R.obs.brn = 0; % temporarily!
        XL = [];
        % Unperturbed
        [~,~,feat_sim{1},~,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
        
        % Now find bursts!!
        R.condname = {'1'};
        [R,BB] = compute_BetaBursts_Simulated(R,xsim_ip{1});
        R.BB.thresh_prctile = 75;% o85; tl 80
        BB = compute_BurstThreshold(R,BB,1,0);
        R.BB.minBBlength = 0.5; %o1 tl 1.5; %  Minimum burst period- cycles
        BB.plot.durlogflag = 0;
        BB = defineBetaEvents(R,BB);
        
        % Now decide the timing of the intervention!
        pulseWid = (0.3.*fsamp);
        pulseAmp = 0.05;
        pulseDelay = (0.*fsamp);
        pulseCycle = (1/50).*fsamp;
        pulseDuty = (0.005*fsamp);
        % If rhythmic Stim
        %         pulseKern = zeros(pulseWid,1);
        %         ps = 1;
        %         while ps <= pulseWid
        %         pulseKern(ps:ps+pulseDuty) = pulseAmp;
        %         ps = ps + floor(pulseCycle+ 5.*randn);
        %         end
        pulseKern = ones(pulseWid,1);
        pU = zeros(size(R.IntP.tvec));
        pulseStart = [];
        for seg = 1:numel(BB.segInds{1})-1
            pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
            pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
            pulse_Phi = BB.Phi{1}(4,pulseInds);
            pulseKern = sin(pulse_Phi+(phaseShift(p))); %
            pU(pulseInds) = pulseKern;
        end
        pU = (0.25.*std(uc{1}(:,1))).*pU; %.*pulseAmp;
        uc_ip{2} =  uc_ip{1};
        uc_ip{2}{1}(:,1) = uc_ip{2}{1}(:,1) + pU'; % Give it a cortical pulse
        %         uc_ip{2}{1}(:,1) = pU'; % Give it a cortical pulse
        [~,~,feat_sim{2},~,xsim_ip{2}]  = computeSimData(R,m,uc_ip{2},Pbase,0);
        
        %         figure(100)
        %         a(1) = subplot(2,1,1)
        %         plot(R.IntP.tvec_obs,xsim_ip{1}{1}(1,2:end));
        %         hold on
        %         plot(R.IntP.tvec_obs,xsim_ip{2}{1}(1,2:end));
        %
        %         a(2) = subplot(2,1,2)
        %
        %         plot(R.IntP.tvec_obs,xsim_ip{1}{1}(4,2:end)+5e-6);
        %         hold on
        %         plot(R.IntP.tvec_obs,xsim_ip{2}{1}(4,2:end)+5e-6);
        %         linkaxes(a,'x')
        %         xlim([7 8])
        %         R.plot.outFeatFx({feat_sim{1}},{feat_sim{2}},R.data.feat_xscale,R,1,[])
        
        % Recompute bursts!
        R.condname = {'1','2'};
        [R,BB] = compute_BetaBursts_Simulated(R,{xsim_ip{1}{1} xsim_ip{2}{1}});
        R.BB.thresh_prctile = 75;% o85; tl 80
        BB = compute_BurstThreshold(R,BB,1,0);
        R.BB.minBBlength = 0.5; %o1 tl 1.5; %  Minimum burst period- cycles
        BB.plot.durlogflag = 0;
        BB = defineBetaEvents(R,BB);
        
        for stm = 1:2
            spec = [squeeze(feat_sim{stm}(1,1,1,1,:)) squeeze(feat_sim{stm}(1,4,4,1,:))];
            powspec_save(:,:,stm,p) = spec;
            intpow(:,1,stm,p) = sum(spec(R.frqz>14 & R.frqz<=21,:));
            maxpow(:,1,stm,p) = max(spec(R.frqz>14 & R.frqz<=21,:));
            intpow(:,2,stm,p) = sum(spec(R.frqz>21 & R.frqz<=30,:));
            maxpow(:,2,stm,p) = max(spec(R.frqz>21 & R.frqz<=30,:));
            
            burdur(stm,p) = median(BB.segDur{stm});
            burAmp(stm,p) = median(BB.segAmp{stm});
            burPPC(stm,p) = median(BB.segPPC{stm});
        end
        
        % Compute Amplitude Corrs etc
        
        % Recompute the burst Envelopes
        null = nan(1,numel(BB.segInds{1})-1);
        XH_save = []; YH_save = []; corrAmp = null; corrPhi = null; xcorrAmp = null; xcorrLAmp = null; nmis = null;
        for seg = 1:numel(BB.segInds{1})-1
            pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
            pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
            if pulseInds(end)<=size(BB.Phi{1},2)
                X = xsim_ip{2}{1}(1,pulseInds);
                Y = xsim_ip{2}{1}(4,pulseInds);
                
                XH = hilbert(X); YH = hilbert(Y);
                XH_save(:,seg) = abs(XH);
                YH_save(:,seg) = abs(YH);
                ctmp = corrcoef(abs(XH),abs(YH));
                corrAmp(seg) = ctmp(1,2);
                RP = angle(XH)-angle(YH);
                corrPhi(seg) =  abs(mean(exp(-1i*RP'),1));
                [c lag] = xcorr(abs(XH),abs(YH),'coeff');
                xcorrAmp(seg) = max(c);
                xcorrLAmp(seg) = lag(c==max(c));
                mid = median([XH YH]);
                %                     nmis(seg) = nmi(XH<mid,YH<mid);
            end
        end
        %             XCor(cond,p) = mean(xcorrAmp);
        %             AEC(cond,p) = mean(corrAmp);
        %             ctmp = corrcoef(mean(XH_save,2),mean(YH_save,2));
        %             AEC_mean(cond,p) = ctmp(1,2);
        %             PLV(cond,p) = mean(corrPhi);
        
        XCor(p) = mean(xcorrAmp);
        AEC(p) = mean(corrAmp);
        ctmp = corrcoef(mean(XH_save,2),mean(YH_save,2));
        AEC_mean(p) = ctmp(1,2);
        PLV(p) = mean(corrPhi);
        
    end
    save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\ClosedLoop_save.mat'],'AEC','AEC_mean','PLV','intpow','powspec_save','maxpow','ck_1','phaseShift')
    
else
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\ClosedLoop_save.mat'],'AEC','AEC_mean','PLV','intpow','powspec_save','maxpow','ck_1','phaseShift')
end

%% Now Plot Results
cmap = brewermap(numel(phaseShift),'Reds');
subplot(1,3,1) % Spectra Plot
phsel = 1:3:18;
ip = 0;
a = [];
for i = phsel
    ip = ip+1;
    a(ip) = plot(R.frqz,1e7.*squeeze(powspec_save(:,2,2,i)),'color',cmap(i,:),'LineWidth',2);
    hold on
end
plot(R.frqz,1e7.*squeeze(powspec_save(:,2,1,1)),'color',[0 0 0],'LineWidth',2,'LineStyle','--');
xlim([8 34])
legend(a,sprintfc('%.1f rad.',phaseShift(phsel)))
xlabel('Frequency (Hz)'); ylabel('STN Power (uV^2 Hz^{-1})')
title('STN Power under Phasic M2 Stimulation')
grid on

subplot(1,3,2) % Steady State Stats
% Beta 1
X = squeeze(intpow(2,1,:,:))';
X = 100.*(X(:,2)-X(:,1))./X(:,1);
p(1) = plot(phaseShift,X,'color',cmap(9,:),'LineWidth',2);
hold on
s = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');

% Beta 2
X = squeeze(intpow(2,2,:,:))';
X = 100.*(X(:,2)-X(:,1))./X(:,1);
p(2) = plot(phaseShift,X,'color',cmap(16,:),'LineWidth',2);
hold on
s(2) = scatter(phaseShift(phsel),X(phsel),75,cmap(phsel,:),'filled');
s(2).Marker = 'square';
grid on

xlabel('Stimulation Phase (radians)'); ylabel('Percentage Change')
legend(p,{'\beta_1 (14-21 Hz) Power','\beta_2 (21-30 Hz) Power'})
title('Spectral Power under Phasic M2 Stimulation')

subplot(1,3,3) % Burst stats
X = AEC;
p(1) = plot(phaseShift,X,'color',cmap(6,:),'LineWidth',2);
hold on
s(1) = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');

X = PLV;
p(2) = plot(phaseShift,X,'color',cmap(12,:),'LineWidth',2);
hold on
s(2) = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');

xlabel('Stimulation Phase (radians)'); ylabel('Percentage Change')
legend(p,{'Burst AEC','Burst PLV'})
title('Burst Statistics under Phasic M2 Stimulation')
grid on

a = gca;
a.XTick = [0 pi/2 pi 3*pi/2 2*pi];
set(gcf,'Position',[711         604        1081         374])

% % % Burst Duration
% % X = burdur';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(1) = plot(phaseShift,X,'color',cmap(6,:),'LineWidth',2);
% % hold on
% % s(1) = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel,:),'filled');
% % 
% % % Burst Amplitude
% % X = burAmp';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(2) = plot(phaseShift,X,'color',cmap(12,:),'LineWidth',2);
% % hold on
% % s(2) = scatter(phaseShift(phsel),X(phsel),65,cmap(phsel,:),'filled');
% % s(2).Marker = 'diamond';
% % 
% % % Burst Phase Locking
% % X = burPPC';
% % X = 100.*(X(:,2)-X(:,1))./X(:,1);
% % p(3) = plot(phaseShift,X,'color',cmap(18,:),'LineWidth',2);
% % hold on
% % s(3) = scatter(phaseShift(phsel),X(phsel),65,cmap(phsel,:),'filled');
% % s(3).Marker = 'square';
% % 
