function [R,MP] = BAA_sim_betaInputs(R,modID,simtime)

cmap = brewermap(18,'Spectral');
condcmap = cmap([1 4 8 16 4 18],:);

% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
MP.p = permMod{1}.par_rep{1};
MP.m = m;
% uc = innovate_timeseries(R,m);
% uc{1} = uc{1}.*sqrt(R.IntP.dt);
% %
% pki = linspace(-pi,pi,20);
% conStren = [0.1 1 1.15];
% for cond = 1:3
%     Pbase = MP.p;
%     Pbase.A{1}(2,1) = -32;
%     Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*conStren(cond)); %
%     parfor i = 1:numel(pki)
%         uc_i = uc;
%         uc_i{1}(:,1) = (0.1.*std(uc{1}(:,1)).*sin(2.*20.*pi.*R.IntP.tvec)) + uc{1}(:,1)';
%         uc_i{1}(:,2) = (0.1.*std(uc{1}(:,2)).*sin(2.*20.*pi.*R.IntP.tvec + pki(i))) + uc{1}(:,2)';
%         
%         
%         [r2mean,pnew,feat_sim,dum,xsim_gl] = computeSimData(R,m,uc_i,Pbase,0);
%         
%         ts_lp{i,cond} = xsim_gl;
%         feat_lp{i,cond} = feat_sim;
%         
%         stn_spec = squeeze(feat_sim(1,4,4,1,:));
%         stn_powspec(:,i,cond) = stn_spec;
%         stn_intpow(i,cond) = sum(stn_spec(R.frqz>14 & R.frqz<25));
%         stn_maxpow(i,cond) = max(stn_spec(R.frqz>14 & R.frqz<25));
%     end
%     
% end
% stn_maxpow_dmin = (stn_maxpow - min(stn_maxpow)); %./std(stn_maxpow);
% stn_maxpow_nmz = 100.*(stn_maxpow - median(stn_maxpow))./median(stn_maxpow);%./std(stn_maxpow);
% for cond = 1:3
%     subplot(1,2,1)
%     b = plot(R.frqz,squeeze(stn_powspec(:,:,cond))');
%     for u = 1:size(b,1)
%         b(u).Color = condcmap(cond+3,:); %.*(1-(u-1)*0.035);
%         b(u).LineWidth = 1;
%     end
%     hold on
%     xlabel('Frequency (Hz)'); ylabel('STN Power')
%     ylim([4e-16,4e-13]); xlim([6 38])
%     set(gca, 'YScale', 'log')
%     grid on
%     subplot(1,2,2)
%     c = plot(pki,stn_maxpow_nmz(:,cond));
%     c.Color =condcmap(cond+3,:);
%     c.LineWidth = 2;
%     hold on
%     xlabel('Relative Phase (\phi_{M2} - \phi_{STR})'); ylabel('% Change in STN Beta Power')
%     xlim([-pi pi])
%     grid on
% %     set(gca, 'YScale', 'log')
% end
% 
% set(gcf,'Position',[651 367 1136 464])
%  R = simannealsetup_InDirect_ModelComp();
% mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
% save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat.mat'],'feat_HD','feat_STR_GPe')
% save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')