function [Rorg] = BAA_sim_BetaPropagation(Rorg,simtime,fresh)
close all
cmap = brewermap(24,'Spectral');
PRC.condcmap = cmap; %([1 4 8 16 4 18],:);
PRC.condname = {'Fitted','1% STN->GPe','150% STN->GPe'}; %Fitted','1% M2->STN','150% M2->STN',
coupname = {'M2/Thal.','STR/Thal.','STN/Thal.'};
N = 12;
% connection list
fsamp = 1/Rorg.IntP.dt;
% ck_1 = logspace(-0.75,0.75,18);
phaseShift = linspace(-pi,pi,N);
for connection = 2
    if connection ==1
        ck_1 = linspace(0.005,3.5,24); % [0.1 1 1.15]; STN-> GPe
        %         ck_1 = logspace(log10(0.05),log10(3.5),18);
    elseif connection == 2
        ck_1 = linspace(0.005,1.5,24); % [0.1 1 1.15]; M2 -> STN
        %         ck_1 = logspace(log10(0.05),log10(1.5),18);
    end
    if fresh
        % Comopute simulations by sweeping across data
        [Rorg,m,permMod,xsimMod{1}] = getSimModelData_v2(Rorg,10,simtime);
        P = permMod{1}.par_rep{1};
        % Give all timeseries the same input - makes comparable
        Rorg.obs.csd.df = 0.25;
        Rorg = setSimTime(Rorg,simtime);
        Rorg.obs.brn = 3; % temporarily!
        uc = innovate_timeseries(Rorg,m);
        uc{1} = uc{1}.*sqrt(Rorg.IntP.dt);
        
        null = cell(numel(ck_1),numel(phaseShift));
        intpow = null;
        maxpow = null;
        XCor = null;
        AEC = null;
        AEC_mean = null;
        PLV = null;
        TE_mean = null;
        ZTE_mean = null;
        Pv_mean = null;
        anTE_mean = null;
        peakTau_mean = null;
        for cond = 1:numel(ck_1)
            parfor p = 1:numel(phaseShift) %[1 10] %
                uc_ip = {};
                uc_ip{1} = uc;
                uc_ip{2} = uc;
                % Setup pulses for PRC computation
                R = Rorg;
                R.obs.brn = 0;
                % Setup the simulations
                Pbase = P;
                
                if connection == 1
                    Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(cond)); %
                else
                    Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(cond)); %
                end
                
                % Unperturbed
                xsim_ip = {}; feat_sim = {};
                [~,~,feat_sim{1},~,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
                
                % Now find bursts!!
                R.condname = {'1'};
                [R,BB] = compute_BetaBursts_Simulated(R,xsim_ip{1});
                R.BB.thresh_prctile = 75;% o85; tl 80
                BB = compute_BurstThreshold(R,BB,1,0);
                R.BB.minBBlength = 2; %o1 tl 1.5; %  Minimum burst period- cycles
                BB.plot.durlogflag = 0;
                BB = defineBetaEvents(R,BB);
                
                % Now decide the timing of the intervention!
                pulseWid = (0.3.*fsamp);
                pulseDelay = (0.*fsamp);
                % If rhythmic Stim
                %         pulseKern = zeros(pulseWid,1);
                %         ps = 1;
                %         while ps <= pulseWid
                %         pulseKern(ps:ps+pulseDuty) = pulseAmp;
                %         ps = ps + floor(pulseCycle+ 5.*randn);
                %         end
                pulseKern = ones(pulseWid,1);
                %             pulseAmp = [-300:0 -1:-1:-299];
                %             pulseAmp = (pulseAmp-min(pulseAmp))./(max(pulseAmp)-min(pulseAmp))
                pU = zeros(size(R.IntP.tvec));
                pulseStart = [];
                for seg = 1:numel(BB.segInds{1})-1
                    pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
                    pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
                    if pulseInds(end)<=size(BB.Phi{1},2)
                        pulse_Phi = BB.Phi{1}(4,pulseInds); %shifted relative to STN phase
                        pulseKern = sin(pulse_Phi+(phaseShift(p))); % shifted
                        pU(pulseInds) = pulseKern;
                    end
                end
                pU = (0.5.*std(uc{1}(:,1))).*pU; %.*pulseAmp;
                uc_ip{2} =  uc_ip{1};
                uc_ip{2}{1}(:,1) = uc_ip{2}{1}(:,1) + pU'; % Give it a cortical pulse
                [~,~,feat_sim{2},~,xsim_ip{2}]  = computeSimData(R,m,uc_ip{2},Pbase,0);
                
                XAct = xsim_ip{2}{1};
                %                 XAct = ft_preproc_bandstopfilter(XAct,2000,[14 34]);
                %                 XActLP = ft_preproc_lowpassfilter(XAct,2000,14);
                %                 XActHB = ft_preproc_highpassfilter(XAct,2000,34);
                
                % Recompute the burst Envelopes
                XH_save = []; YH_save = [];
                corrAmp_0 = []; corrPhi_0 = []; xcorrAmp_0= []; xcorrLAmp_0 = []; nmis = [];
                TE_0 = []; Pv_0 = []; anTE_0 = []; peakTau_0 = []; ZTE_0 = [];
                corrAmp_1 = []; corrPhi_1 = []; xcorrAmp_1= []; xcorrLAmp_1 = []; nmis = [];
                TE_1 = []; Pv_1 = []; anTE_1 = []; peakTau_1 = []; ZTE_1 = [];
                for seg = 1:numel(BB.segInds{1})-1
                    pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
                    pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
                    if pulseInds(end)<=size(BB.Phi{1},2)
                        X1 = xsim_ip{1}{1}(1,pulseInds);
                        X2 = xsim_ip{1}{1}(2,pulseInds);
                        X3 = xsim_ip{1}{1}(4,pulseInds);
                        X = reshape([X1;X2;X3],3,1,[]);
                        Y = xsim_ip{1}{1}(5,pulseInds);
                        
                        [corrAmp_0(seg,:,:),xcorrAmp_0(seg,:,:),corrPhi_0(seg,:,:),xcorrLAmp_0(seg,:,:),TE_0(seg,:,:,:),...
                            Pv_0(seg,:,:,:),anTE_0(seg,:,:),peakTau_0(seg,:,:,:),ZTE_0(seg,:,:,:)] = computePropMetrics(X,Y);
                        
                        X1 = xsim_ip{2}{1}(1,pulseInds);
                        X2 = xsim_ip{2}{1}(2,pulseInds);
                        X3 = xsim_ip{2}{1}(3,pulseInds);
                        X = reshape([X1;X2;X3],3,1,[]);
                        Y = xsim_ip{2}{1}(4,pulseInds);
                        
                        [corrAmp_1(seg,:,:),xcorrAmp_1(seg,:,:),corrPhi_1(seg,:,:),xcorrLAmp_1(seg,:,:),TE_1(seg,:,:,:),...
                            Pv_1(seg,:,:,:),anTE_1(seg,:,:),peakTau_1(seg,:,:,:),ZTE_1(seg,:,:,:)] = computePropMetrics(X,Y);
                        
                        
                        
                        
                        %                         X1 = XActLP(1,pulseInds);
                        %                         X2 = XActLP(2,pulseInds);
                        %                         X3 = XActLP(4,pulseInds);
                        %                         X4 = XActHB(1,pulseInds);
                        %                         X5 = XActHB(2,pulseInds);
                        %                         X6 = XActHB(4,pulseInds);
                        %                         X = reshape([X1;X2;X3;X4;X5;X6],3,2,[]);
                        %
                        %                         Y1 =  XActLP(5,pulseInds);
                        %                         Y2 =  XActHB(5,pulseInds);
                        %                         Y = [Y1;Y2];
                        
                        % TE(1) = STN -> M2
                        % TE(2) = M2 -> STN
                        %                     nmis(seg) = nmi(XH<mid,YH<mid);
                    end
                end
                XCor{cond,p} = mean(xcorrAmp);
                AEC{cond,p} = mean(corrAmp);
                ctmp = corrcoef(mean(XH_save,2),mean(YH_save,2));
                AEC_mean{cond,p} = ctmp(1,2);
                PLV{cond,p} = mean(corrPhi);
                TE_mean{cond,p} = mean(TE);
                ZTE_mean{cond,p} = mean(ZTE);
                Pv_mean{cond,p} = mean(Pv);
                anTE_mean{cond,p} = mean(anTE);
                peakTau_mean{cond,p} = (1000.*mean(peakTau))./fsamp;
                intpow_tmp = []; maxpow_tmp = []; powspec_tmp = [];
                for stm = 1:2
                    
                    spec = [squeeze(feat_sim{stm}(1,1,1,1,:)) squeeze(feat_sim{stm}(1,4,4,1,:))];
                    powspec_tmp(:,:,stm) = spec;
                    intpow_tmp(:,1,stm) = sum(spec(R.frqz>14 & R.frqz<=21,:));
                    maxpow_tmp(:,1,stm) = max(spec(R.frqz>14 & R.frqz<=21,:));
                    intpow_tmp(:,2,stm) = sum(spec(R.frqz>21 & R.frqz<=30,:));
                    maxpow_tmp(:,2,stm) = max(spec(R.frqz>21 & R.frqz<=30,:));
                end
                
                intpow{cond,p} = intpow_tmp; % channel band stim K phi
                maxpow{cond,p} = maxpow_tmp;
                powspec{cond,p} = powspec_tmp;
            end
        end
        
        
        if connection == 1
            save([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_HD.mat'],'XCor','AEC','PLV','intpow','powspec','maxpow','ck_1','phaseShift','TE_mean','Pv_mean','anTE_mean','peakTau_mean','ZTE_mean')
        elseif connection == 2
            save([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_STNGPe.mat'],'XCor','AEC','PLV','intpow','powspec','maxpow','ck_1','phaseShift','TE_mean','Pv_mean','anTE_mean','peakTau_mean','ZTE_mean')
        end
    else
        if connection == 1
            load([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_HD.mat'],'XCor','AEC','PLV','intpow','powspec','maxpow','ck_1','phaseShift','TE_mean','Pv_mean','anTE_mean','peakTau_mean','ZTE_mean')
        elseif connection == 2
            load([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_STNGPe.mat'],'XCor','AEC','PLV','intpow','powspec','maxpow','ck_1','phaseShift','TE_mean','Pv_mean','anTE_mean','peakTau_mean','ZTE_mean')
        end
    end
    phaseShift = linspace(-pi,pi,N);
    
    for i = 1:8
        if i == 1
            LP = unPackCell(XCor);
            titname = 'Envelope Correlation';
            cl = [0.8 0.95];
            dcl = [-0.05 0];
        elseif i == 2
            LP = unPackCell(XCor);
            titname = 'Phase Locking';
            cl = [0.25 0.8];
            dcl = [-0.1 0.04];
        elseif i == 3
            LP = unPackCell(TE1);
            titname = 'TE1';
            cl = [0.01 0.04];
            dcl = [-0.05 0];
            %LP(LP==0) = NaN;
            %         LP = LP.*squeeze(Pv_mean(1,:,:,:)<0.05);
        elseif i == 4
            LP = squeeze(TE_mean(2,:,:,:)); titname = 'TE2';
            cl = [0.015 0.04];
            dcl = [-0.01 0.01];
            %         LP = LP.*squeeze(Pv_mean(2,:,:,:)<0.05);
        elseif i == 5
            LP = squeeze(Pv_mean(1,:,:,:)); titname = 'PV'; LP(LP==0) = NaN;
            cl = [0.025 0.04];
            dcl = [-0.01 0.01];
            %             LP = LP.*squeeze(Pv_mean(2,:,:)<0.05);
        elseif i == 6
            LP = anTE_mean; titname = 'deltaTE';
            cl = [-0.35 0.15];
            dcl = [-0.1 0.5];
            %             LP = LP.*squeeze(Pv_mean(2,:,:)<0.1);
            %             LP(LP==0)    = -32;
            %             crng(crng<0.01
        elseif i == 7
            LP = squeeze(ZTE_mean(1,:,:,:)); titname = 'STN -> M2';
            cl = [0.025 0.04];
            dcl = [-0.01 0.01];
            %             LP = LP.*squeeze(Pv_mean(1,:,:)<0.1);
            %             crng = [0 linspace(cl(1),cl(2),9)];
        elseif i == 8
            LP = squeeze(ZTE_mean(2,:,:,:)); titname = 'M2 -> STN';
            cl = [0.025 0.04];
            dcl = [-0.01 0.01];
            %             LP = LP.*squeeze(Pv_mean(2,:,:)<0.1);
            %             crng = [0 linspace(cl(1),cl(2),9)];
        end
        
        crng = linspace(cl(1),cl(2),12);
        dcrng = linspace(dcl(1),dcl(2),12);
        
        LP(isnan(AEC)) = NaN;
        figure((connection*10)+i)
        for p= 1:3
            subplot(1,3,p)
            %  [Fig] = CirHeatmap({LP'}, 'GroupLabels', '1','OuterLabels',sprintfc('%.f',rad2deg(phaseShift)), 'CircType', 'o','InnerSpacerSize',0,'LogRadius',1);
            %  caxis(cl)
            %         newpoints = 100;
            %         [xq,yq] = meshgrid(...
            %             linspace(min(min(phaseShift,[],2)),max(max(phaseShift,[],2)),newpoints ),...
            %             linspace(min(min(log10(ck_1),[],1)),max(max(log10(ck_1),[],1)),newpoints )...
            %             );
            %         LPint = interp2(phaseShift,log10(ck_1),LP,xq,yq,'cubic');
            %         [a b c] = contourf(xq,yq,LPint,crng)
            if p<3
                contourf(phaseShift,ck_1*100,squeeze(LP(p,:,:)),crng);
                Q = cl;
            else
                contourf(phaseShift,ck_1*100,squeeze(LP(2,:,:))-squeeze(LP(1,:,:)),dcrng);
                Q = dcl;
            end
        end
        
        % Compute Power Responses
        stn_maxpow = squeeze(maxpow(2,1,2,:,:))';
        stn_maxpow_dmin = (stn_maxpow - min(stn_maxpow)); %./std(stn_maxpow);
        stn_maxpow_nmz = 100.*(stn_maxpow - median(stn_maxpow))./median(stn_maxpow);%./std(stn_maxpow);
        condcmap = brewermap(24,'Reds');
    end
    
end

%% SCRIPT GRAVE
%             %  [Fig] = CirHeatmap({LP'}, 'GroupLabels', '1','OuterLabels',sprintfc('%.f',rad2deg(phaseShift)), 'CircType', 'o','InnerSpacerSize',0,'LogRadius',1);
%             %  caxis(cl)
%             %         newpoints = 100;
%             %         [xq,yq] = meshgrid(...
%             %             linspace(min(min(phaseShift,[],2)),max(max(phaseShift,[],2)),newpoints ),...
%             %             linspace(min(min(log10(ck_1),[],1)),max(max(log10(ck_1),[],1)),newpoints )...
%             %             );
%             %         LPint = interp2(phaseShift,log10(ck_1),LP,xq,yq,'cubic');
%             %         [a b c] = contourf(xq,yq,LPint,crng)

% surf(phaseShift,log10(ck_1),squeeze(PLV(3,:,:))')
%     figure
%     for cond = 1:2:18
%         subplot(1,3,1)
%         b = plot(Rorg.frqz,squeeze(powspec(:,2,2,cond,:))');
%         hold on
%         for u = 1:size(b,1)
%             b(u).Color = condcmap(cond+6,:); %.*(1-(u-1)*0.035);
%             b(u).LineWidth = 1;
%         end
%         hold on
%         xlabel('Frequency (Hz)'); ylabel('STN Power')
%         ylim([4e-16,4e-13]); xlim([6 38])
%         set(gca, 'YScale', 'log')
%         grid on
%
%         subplot(1,3,2)
%         c = plot(phaseShift,stn_maxpow_nmz(:,cond));
%         c.Color =condcmap(cond+6,:);
%         c.LineWidth = 2;
%
%         hold on
%         xlabel('Relative Phase (\phi_{M2} - \phi_{STR})'); ylabel('% Change in STN Beta Power')
%         xlim([0 2*pi])
%         grid on
%     end
%
%
%     subplot(1,3,3)
%     p = plotPRCSumStats(ck_1,max(stn_maxpow_nmz),min(stn_maxpow_nmz),range(stn_maxpow_nmz),1:2:18,condcmap);
%
%     set(p,'Color',cmap(end,:),'LineWidth',2)
%     ylabel('% Change in STN Beta Amplitude')
%     xlabel('Connection Strength (% of fitted)')
%     grid on
%     legend(p,'PRC Max','PRC Min','PRC Range')
%     set(gca, 'XScale', 'log')
%     set(gcf,'Position',[103 595 1362 360])
