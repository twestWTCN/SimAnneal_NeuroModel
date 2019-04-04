close all
% rng(1231)
% RESULTS ARENT REPRODUCIBLE AT ALL - IS THIS LARGE AMOUNTS OF NOISE? CUT
% OUT SOME NOISE!
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
i = randi(25);
varo = load([R.rootn '\Inputs\TMSsig\TMSsignal_' num2str(i)]);
u = varo.varo;
j = randi(25);
varo = load([R.rootn '\Inputs\betaSig\betaSignal_' num2str(j)]);
m = varo.varo;

% I2amp = [-2 0 2];
% betaSNR = [-1 0 1]; % The background beta power decides the slope of the
% amp vs latency (sensitive parameter!)
CSNSNR = [-2 0 2];
 
for isweep = 1:3
    
    pfit = spm_unvec(parBank(1:end-1,1),pt);
%     pfit.IWS_amp = [1 0];
%     pfit.IWS_amp_jit = [0 -2];
%     pfit.SNRs = [0.2 CSNSNR(isweep) 0.4 -0.5];
%     pfit.EPSP_amp = [0.25 1];
% p.CSN2AMN = CSNSNR(isweep);
    % [r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,BPfit,0,1);
    R = PB_scheme_setup();
    [R] = getFlaviesData(R);
    % pfit.EPSP_ampJit = [0 0 0];
    [xsims,dum,wflag] = PB_schema_simulate_REV(R,[],u,pfit,m);
    
    [dum,data_emp] = getFlaviesData(R);
    
    [~, feat_sim, wflag(2)] = R.obs.transFx(R,xsims);
    figure(1);  R.plot.outFeatFx(R,{R.data.feat_emp},{feat_sim},1)
    
    amp_emp = data_emp.AmpvsLat(:,1);
    lat_emp = data_emp.AmpvsLat(:,2)*1000;
    amp_sim = xsims.zMEP_amp;
    lat_sim = xsims.MEP_mconset;
    % lat_sim = xsims.zMEP_onset;
    lat_supp = -3:.25:3;
    
    figure(2)
    subplot(2,1,1)
    histogram(lat_sim,lat_supp,'Normalization','probability')
    hold on
    %   histogram(lat_emp,lat_supp,'Normalization','probability')
    xlabel('Mean Corrected Latency (ms)'); title('All MEP')
    subplot(2,2,3)
    histogram(lat_sim(amp_sim>1.2),lat_supp,'Normalization','probability')
    xlabel('Mean Corrected Latency (ms)'); title('MEP > 1.2z')
      hold on
    %    histogram(lat_emp(amp_emp>1.2),lat_supp,'Normalization','probability')
    subplot(2,2,4)
    histogram(lat_sim(amp_sim<1.2),lat_supp,'Normalization','probability')
    hold on
    %    histogram(lat_emp(amp_emp<1.2),lat_supp,'Normalization','probability')
    xlabel('Mean Corrected Latency (ms)'); title('MEP < 1.2z')
    
    figure(3)
    scatter(amp_sim,lat_sim);
    hold on
    % scatter(amp_emp,lat_emp)
    xlabel('MEP Amplitude (z-scored)'); ylabel('Mean Corrected Latency')
    
    figure(4)
    subplot(2,1,1)
    scatter(xsims.TMS_phase,amp_sim)
    hold on
    [binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,amp_sim);
    binStats_phizAmp(end,:) = [];
    plot(binEdge2Mid(R.trans.pirange),binStats_phizAmp(:,3)')
    xlabel('Pre-TMS Beta Phase'); ylabel('MEP Amplitude (z-scored)');
    subplot(2,1,2)
    scatter(xsims.TMS_phase,lat_sim)
    hold on
    [binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,lat_sim);
    binStats_phizAmp(end,:) = [];
    plot(binEdge2Mid(R.trans.pirange),binStats_phizAmp(:,3)')
    xlabel('Pre-TMS Beta Phase'); ylabel('Mean Corrected Latency (ms)');
    figure(1)
end
