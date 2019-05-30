function [R] = BAA_sim_betaInputs(R,modID,simtime)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
XBase = permMod{1}.par_rep{1};

% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);

pki = linspace(-pi,pi,10);
parfor i = 1:numel(pki)
    uc_i = uc;
    uc_i{1}(:,1) = (0.1.*std(uc{1}(:,1)).*sin(2.*20.*pi.*R.IntP.tvec)) + uc{1}(:,1)';
    uc_i{1}(:,2) = (0.1.*std(uc{1}(:,2)).*sin(2.*20.*pi.*R.IntP.tvec + pki(i))) + uc{1}(:,2)';
    
    Pbase = XBase;
    Pbase.A{1}(2,1) = -32;
    [r2mean,pnew,feat_sim,dum,xsim_gl] = computeSimData(R,m,uc_i,Pbase,0);
    
    ts_lp{i} = xsim_gl;
    feat_lp{i} = feat_sim;
    
    stn_spec = squeeze(feat_sim(1,4,4,1,:));
    stn_powspec(:,i) = stn_spec;
    stn_intpow(i) = sum(stn_spec(R.frqz>14 & R.frqz<25));
    stn_maxpow(i) = max(stn_spec(R.frqz>14 & R.frqz<25));
end
subplot(1,2,1)
plot(R.frqz,stn_powspec')
xlabel('Frequency (Hz)'); ylabel('STN Power')
xlim([4 3   2])
subplot(1,2,2)
plot(pki,stn_maxpow)
xlabel('Relative Phase (\phi_{M2} - \phi_{STR})'); ylabel('STN Power')
xlim([-pi pi])
R = simannealsetup_InDirect_ModelComp();
mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat.mat'],'feat_HD','feat_STR_GPe')
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')