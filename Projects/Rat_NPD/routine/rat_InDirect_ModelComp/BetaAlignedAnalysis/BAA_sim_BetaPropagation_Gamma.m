function [Rorg] = BAA_sim_BetaPropagation_Gamma(Rorg,simtime,fresh)
close all
cmap = brewermap(24,'Spectral');
PRC.condcmap = cmap; %([1 4 8 16 4 18],:);
PRC.condname = {'Fitted','1% STN->GPe','150% STN->GPe'}; %Fitted','1% M2->STN','150% M2->STN',
coupname = {'M2/Thal.','STR/Thal.','STN/Thal.'};
N = 6;
% connection list
fsamp = 1/Rorg.IntP.dt;
% ck_1 = logspace(-0.75,0.75,18);
phaseShift = linspace(-pi,pi,N);
for connection = 1:2
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
        [t,gmx,gmxAmp,gmxPhi] = makeBurstSignal(Rorg.IntP.tvec,70,3);
        gmx = 5.*gmx;
        uc = innovate_timeseries(Rorg,m);
        uc{1} = uc{1}.*sqrt(Rorg.IntP.dt);
        intpow = nan(2,2,2,numel(ck_1),numel(phaseShift)); maxpow = nan(2,2,2,numel(ck_1),numel(phaseShift));
        
        for cond = 12; %1:numel(ck_1)
            for p = 1:numel(phaseShift) %[1 10] %
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
                uc_ip{1}{1}(2:end,1) = uc_ip{1}{1}(2:end,1) + gmx';
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
                        pU(pulseInds) = pulseKern ;
                    end
                end
                pU = (0.5.*std(uc{1}(:,1))).*pU; %.*pulseAmp;
                uc_ip{2} =  uc_ip{1};
                uc_ip{2}{1}(:,1) = uc_ip{2}{1}(:,1) + pU'; % Give it a cortical pulse
                [~,~,feat_sim{2},~,xsim_ip{2}]  = computeSimData(R,m,uc_ip{2},Pbase,0);
                
                XAct = xsim_ip{1}{1};
                XActHB_0 = ft_preproc_bandpassfilter(XAct,2000,[60 80]);
                XAct = xsim_ip{2}{1};
                XActHB_1 = ft_preproc_bandpassfilter(XAct,2000,[60 80]);                
                %                 XAct = ft_preproc_bandstopfilter(XAct,2000,[14 34]);
%                 XActLP = ft_preproc_lowpassfilter(XAct,2000,14);
%                 XActHB = ft_preproc_highpassfilter(XAct,2000,40);
                % Recompute the burst Envelopes
                XH_save = []; YH_save = []; xcorrAmp_0 = []; xcorrAmp_1 = []; corrPhi_0 = []; corrPhi_1 = []; xcorrAmp = []; xcorrLAmp = []; nmis = [];
                TE = []; Pv = []; anTE = []; peakTau = []; ZTE = [];
                for seg = 1:numel(BB.segInds{1})-1
                    pulseStart(seg) = BB.segInds{1}{seg}(1) + pulseDelay;
                    pulseInds = pulseStart(seg):pulseStart(seg)+pulseWid-1;
                    if pulseInds(end)<=size(BB.Phi{1},2)
                        X1 = XActHB_0(1,pulseInds);
                        X2 = XActHB_0(2,pulseInds);
                        X3 = XActHB_0(3,pulseInds);
                        X = reshape([X1;X2;X3],3,1,[]);
                        Y  = XActHB_0(4,pulseInds);
                        [dum,xcorrAmp_0(:,seg),corrPhi_0(:,seg)] = computePropMetrics(X,Y);
                        
                        X1 = XActHB_1(1,pulseInds);
                        X2 = XActHB_1(2,pulseInds);
                        X3 = XActHB_1(3,pulseInds);
                        X = reshape([X1;X2;X3],3,1,[]);
                        Y  = XActHB_1(4,pulseInds);
                        [dum,xcorrAmp_1(:,seg),corrPhi_1(:,seg)] = computePropMetrics(X,Y);
                        %                     nmis(seg) = nmi(XH<mid,YH<mid);
                    end
                end
                XCor(:,cond,p) = squeeze(mean(xcorrAmp_1'));
                dXCor(:,cond,p) = squeeze(mean((xcorrAmp_0-xcorrAmp_1)'));
                PLV(:,cond,p) = squeeze(mean(corrPhi_1'));
                dPLV(:,cond,p) = squeeze(mean((corrPhi_0-corrPhi_1)'));
                intpow_tmp = []; maxpow_tmp = []; powspec_tmp = [];
                for stm = 1:2
                    spec = [squeeze(feat_sim{stm}(1,1,1,1,:)) squeeze(feat_sim{stm}(1,4,4,1,:))];
                    powspec_tmp(:,:,stm) = spec;
                    intpow_tmp(:,1,stm) = sum(spec(R.frqz>14 & R.frqz<=21,:));
                    maxpow_tmp(:,1,stm) = max(spec(R.frqz>14 & R.frqz<=21,:));
                    intpow_tmp(:,2,stm) = sum(spec(R.frqz>21 & R.frqz<=30,:));
                    maxpow_tmp(:,2,stm) = max(spec(R.frqz>21 & R.frqz<=30,:));
                end
                
                intpow(:,:,:,cond,p) = intpow_tmp; % channel band stim K phi
                maxpow(:,:,:,cond,p) = maxpow_tmp;
                powspec(:,:,:,cond,p)= powspec_tmp;
                disp([cond p])
            end
        end
        if connection == 1
            save([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_Gamma_HD.mat'],'XCor','dXCor','PLV','dPLV','intpow','powspec','maxpow','ck_1','phaseShift')
        elseif connection == 2
            save([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_Gamma_STNGPe.mat'],'XCor','dXCor','PLV','dPLV','intpow','powspec','maxpow','ck_1','phaseShift')
        end
    else
        if connection == 1
            load([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_Gamma_HD.mat'],'XCor','dXCor','PLV','dPLV','intpow','powspec','maxpow','ck_1','phaseShift')
        elseif connection == 2
            load([Rorg.rootn '\routine\' Rorg.out.tag '\BetaBurstAnalysis\Data\BetaPropagation_Gamma_STNGPe.mat'],'XCor','dXCor','PLV','dPLV','intpow','powspec','maxpow','ck_1','phaseShift')
        end
    end
    phaseShift = linspace(-pi,pi,N);
    
    for bnd = 1:2
        for i = 1:4
            if i == 1
                LP = squeeze(XCor(:,:,:)); titname = 'Envelope Correlation';
                cl = [0.8 0.95];
                dcl = [-0.05 0];
            elseif i == 2
                LP = squeeze(PLV(:,:,:)); titname = 'Phase Locking';
                cl = [0.25 0.8];
                dcl = [-0.1 0.04];
            elseif i == 3
                LP = squeeze(TE_mean(bnd,1,:,:,:)); titname = 'TE For';
                cl = [0.01 0.06];
                dcl = [-0.05 0];
                %LP(LP==0) = NaN;
                %         LP = LP.*squeeze(Pv_mean(1,:,:,:)<0.05);
            elseif i == 4
                LP = squeeze(TE_mean(bnd,2,:,:,:)); titname ='TE Rev';
                cl = [0.015 0.06];
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
                LP = squeeze(ZTE_mean(1,:,:,:)); titname = 'ZTE For';
                cl = [0.025 0.04];
                dcl = [-0.01 0.01];
                %             LP = LP.*squeeze(Pv_mean(1,:,:)<0.1);
                %             crng = [0 linspace(cl(1),cl(2),9)];
            elseif i == 8
                LP = squeeze(ZTE_mean(2,:,:,:)); titname = 'ZTE Rev';
                cl = [0.025 0.04];
                dcl = [-0.01 0.01];
                %             LP = LP.*squeeze(Pv_mean(2,:,:)<0.1);
                %             crng = [0 linspace(cl(1),cl(2),9)];
            elseif i == 9
                LP(1,:,:) = squeeze(intpow(1,1,2,:,:))./squeeze(intpow(1,1,1,:,:)); titlname{1} = 'M2 LB';
                LP(2,:,:) = squeeze(intpow(2,1,2,:,:))./squeeze(intpow(2,1,1,:,:)); titlname{2} = 'STN LB';
                LP(3,:,:) = squeeze(intpow(1,2,2,:,:))./squeeze(intpow(1,2,1,:,:)); titlname{3} = 'M2 HB';
                LP(4,:,:) = squeeze(intpow(2,2,2,:,:))./squeeze(intpow(2,2,1,:,:)); titlname{4} = 'STN HB';
                cl = [70 195];
                dcl = [70 195];
            end
            
            crng = linspace(cl(1),cl(2),12);
            dcrng = linspace(dcl(1),dcl(2),12);
            
%             LP(isnan(AEC)) = NaN;
            figure((connection*10)+i)
            
            if i<9
                for p= 1:3
                    subplot(2,3,p+(bnd-1)*3)
                    if p<4
                        contourf(phaseShift,(ck_1*100),squeeze(LP(p,:,:))); %,crng);
                        Q = cl;
                    else
                        contourf(phaseShift,(ck_1*100),squeeze(LP(2,:,:))-squeeze(LP(1,:,:))) %,dcrng);
                        Q = dcl;
                    end
                    colorbar
                    ylabel('log Coupling Strength'); xlabel('angle'); title([coupname{p} ' ' titname])
                    colormap(flipud(brewermap(128,'RdYlBu')))
                    a = gca;
                    a.XTick = [-pi -pi/2 0 pi/2 pi ];
                    %                 caxis(Q)
                    axis square
                    set(gcf,'Position',[680         541        1125         437])
                    %  a.XTickLabel = {'- \pi','-\pi/2','0','+\pi/2','+\pi'}
                end
            elseif i==9
                for p= 1:4
                    subplot(2,2,p)
                    P = fillmissing(squeeze(LP(p,:,:)),'spline');
                    contourf(phaseShift,(ck_1*100),P*100) %,crng);
                    Q = cl;
                    
                    colorbar
                    ylabel('log Coupling Strength'); xlabel('angle'); title([titlname{p}])
                    colormap(flipud(brewermap(128,'RdYlBu')))
                    a = gca;
                    a.XTick = [-pi -pi/2 0 pi/2 pi ];
                    %                 caxis(Q)
                    axis square
                    
                end
                set(gcf,'Position',[680         541        1125         437])
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
