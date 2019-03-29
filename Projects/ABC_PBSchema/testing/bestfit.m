close all
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
Mfit = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;
[dum pt mt] = setup_model_testing(R);

clear varo
R = PB_scheme_setup();
% precomputeBetaSig(R,mt)
% precomputeTMS(R,mt)
varo = load([R.rootn '\Inputs\TMSsig\TMSsignal_' num2str(randi(25))]);
u = varo.varo;
varo = load([R.rootn '\Inputs\betaSig\betaSignal_' num2str(randi(25))]);
m = varo.varo;

pfit = spm_unvec(parBank(1:end-1,1),pt);
% [r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,BPfit,0,1);
R = PB_scheme_setup();
[R] = getFlaviesData(R);
[xsims,dum,wflag] = PB_schema_simulate_REV(R,[],u,pfit,m);

[dum,data_emp] = getFlaviesData(R);

[~, feat_sim, wflag(2)] = R.obs.transFx(R,xsims);
figure(1);  R.plot.outFeatFx(R,{R.data.feat_emp},{feat_sim},1)

amp_emp = data_emp.AmpvsLat(:,1);
lat_emp = data_emp.AmpvsLat(:,2)*1000;
amp_sim = xsims.zMEP_amp;
%   lat_sim = xsims.MEP_mconset;
lat_sim = xsims.zMEP_onset;
lat_supp = -3:.2:3;

figure
subplot(2,1,1)
histogram(lat_sim,lat_supp,'Normalization','probability')
%   hold on
%   histogram(lat_emp,lat_supp,'Normalization','probability')
% xlabel('Z
subplot(2,2,3)
histogram(lat_sim(amp_sim>1.2),lat_supp,'Normalization','probability')
%   hold on
%    histogram(lat_emp(amp_emp>1.2),lat_supp,'Normalization','probability')
subplot(2,2,4)
histogram(lat_sim(amp_sim<1.2),lat_supp,'Normalization','probability')
%   hold on
%    histogram(lat_emp(amp_emp<1.2),lat_supp,'Normalization','probability')

figure
scatter(amp_sim,lat_sim);
hold on
scatter(amp_emp,lat_emp)
xlabel('MEP Amplitude'); ylabel('MEP Latency')

figure
subplot(2,1,1)
scatter(xsims.TMS_phase,amp_sim)
hold on
[binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,amp_sim);
binStats_phizAmp(end,:) = [];
plot(binEdge2Mid(R.trans.pirange),binStats_phizAmp(:,3)')

subplot(2,1,2)
scatter(xsims.TMS_phase,lat_sim)
hold on
[binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,lat_sim);
binStats_phizAmp(end,:) = [];
plot(binEdge2Mid(R.trans.pirange),binStats_phizAmp(:,3)')

