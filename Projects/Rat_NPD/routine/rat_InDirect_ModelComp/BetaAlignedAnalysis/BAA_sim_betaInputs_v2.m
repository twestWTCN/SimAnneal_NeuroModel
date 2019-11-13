function [R,MP] = BAA_sim_betaInputs(R,modID,simtime)

cmap = brewermap(24,'Reds');
condcmap = cmap; %cmap([1 4 8 16 4 18],:);

% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
MP.p = permMod{1}.par_rep{1};
MP.m = m;
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
%
pki = linspace(-pi,pi,20);
NcS = 18;
for ctype = 2
    if ctype == 1
        conStren = linspace(0.05,1.5,NcS); % [0.1 1 1.15]; STN-> GPe
    elseif ctype == 2
        conStren = linspace(0.05,2.5,NcS); % [0.1 1 1.15]; M2 -> STN
    end
    stn_powspec = [];
    stn_intpow = [];
    stn_maxpow = [];
    
    for cond = 1:NcS
        Pbase = MP.p;
        if ctype ==1
            
            Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*conStren(cond)); %
        elseif ctype == 2
            Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*conStren(cond)); %
        end
        
        %     Pbase.A{1}(2,1) = -32; % Set M2->STR to zero
        %     Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*conStren(cond)); %
        %     Pbase.A{2}(3,2) = log(exp(Pbase.A{2}(3,2))*conStren(cond)); %
        %     Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*conStren(cond)); %
        
        parfor i = 1:numel(pki)
            uc_i = uc;
            uc_i{1}(:,1) = (0.1.*std(uc{1}(:,1)).*sin(2.*20.*pi.*R.IntP.tvec)) + uc{1}(:,1)'; % M2
            %         uc_i{1}(:,2) = (0.1.*std(uc{1}(:,2)).*sin(2.*20.*pi.*R.IntP.tvec + pki(i))) + uc{1}(:,2)';
            uc_i{1}(:,4) = (0.1.*std(uc{1}(:,2)).*sin(2.*20.*pi.*R.IntP.tvec + pki(i))) + uc{1}(:,2)'; % STN
            
            [r2mean,pnew,feat_sim,dum,xsim_gl] = computeSimData(R,m,uc_i,Pbase,0);
            
            
            X1 = xsim_gl{1}(1,:);
            X2 = xsim_gl{1}(2,:);
            X3 = xsim_gl{1}(3,:);
            X = reshape([X1;X2;X3],3,1,[]);
            Y = xsim_gl{1}(4,:);
            
            [corrAmp(i,cond,:),xcorrAmp(i,cond,:),corrPhi(i,cond,:),xcorrLAmp(i,cond,:),TE(i,cond,:,:),...
                Pv(i,cond,:,:),anTE(i,cond,:),peakTau(i,cond,:,:),ZTE(i,cond,:,:)] = computePropMetrics(X,Y);
            
            ts_lp{i,cond} = xsim_gl;
            feat_lp{i,cond} = feat_sim;
            
            stn_spec = squeeze(feat_sim(1,4,4,1,:));
            stn_powspec(:,i,cond) = stn_spec;
            stn_intpow(i,cond) = sum(stn_spec(R.frqz>14 & R.frqz<30));
            stn_maxpow(i,cond) = max(stn_spec(R.frqz>14 & R.frqz<30));
        end
        
    end
    if ctype == 2
        save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BetaM2STRInput_HD.mat'],'corrAmp','xcorrAmp','corrPhi','xcorrLAmp','TE','Pv','anTE','peakTau','ZTE','stn_intpow','stn_maxpow')
    elseif ctype == 1
        save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BetaM2STRInput_STNGPe.mat'],'corrAmp','xcorrAmp','corrPhi','xcorrLAmp','TE','Pv','anTE','peakTau','ZTE','stn_intpow','stn_maxpow')
    end
end
% 
% % load('PRC_betainput_tmp.mat')
% %     stn_maxpow_dmin = (stn_maxpow - min(stn_maxpow)); %./std(stn_maxpow);
% %     stn_maxpow_nmz = 100.*(stn_maxpow - median(stn_maxpow))./median(stn_maxpow);%./std(stn_maxpow);
% %     if ctype ==1
% %         cmap = brewermap(24,'Reds');
% %     elseif ctype == 2
% %         cmap = brewermap(24,'Blues');
% %     end
% %     figure(ctype)
% %     condcmap = cmap;
% %     NcS_sel = [1 5 10 14 18];
% %     for cond = NcS_sel
% %         subplot(1,3,1)
% %         b = plot(R.frqz,squeeze(stn_powspec(:,:,cond))');
% %         for u = 1:size(b,1)
% %             b(u).Color = condcmap(cond+6,:); %.*(1-(u-1)*0.035);
% %             b(u).LineWidth = 1;
% %         end
% %         hold on
% %         xlabel('Frequency (Hz)'); ylabel('STN Power')
% %         ylim([1e-16,4e-13]); xlim([6 38])
% %         set(gca, 'YScale', 'log')
% %         grid on
% %         subplot(1,3,2)
% %         c = plot(pki,stn_maxpow_nmz(:,cond));
% %         c.Color =condcmap(cond+6,:);
% %         c.LineWidth = 2;
% %         a = gca;
% %         a.XTick = [-pi -pi/2 0 pi/2 pi];
% %
% %         hold on
% %         xlabel('Relative Phase (\phi_{M2} - \phi_{STR})'); ylabel('% Change in STN Beta Power')
% %         xlim([-pi pi])
% %         grid on
% %         %     set(gca, 'YScale', 'log')
% %
% %         %     [xq yq R2(cond) s(:,cond)] = sinfit(pki',stn_maxpow_nmz(:,cond),50,0);
% %     end
% %
% %     subplot(1,3,3)
% %     p = plotPRCSumStats(conStren,max(stn_maxpow_nmz),min(stn_maxpow_nmz),range(stn_maxpow_nmz),NcS_sel,condcmap);
% %     set(p,'Color',cmap(end,:),'LineWidth',2)
% %     ylabel('% Change in STN Beta Amplitude')
% %     xlabel('Connection Strength (% of fitted)')
% %     grid on
% %     legend(p,'PRC Max','PRC Min','PRC Range')
% %     xlim([5 150])
% %
% %     set(gcf,'Position',[103 595 1362 360])
% %
% % end


%  R = simannealsetup_InDirect_ModelComp();
% mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
% save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat.mat'],'feat_HD','feat_STR_GPe')
% save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')