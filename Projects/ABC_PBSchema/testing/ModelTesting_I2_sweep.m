close all
cmap = lines(7);
% rng(1231)
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
CSNSNR = [-3 0 3];

for isweep = 1:3
    
    pfit = spm_unvec(parBank(1:end-1,1),pt);
    %         pfit.IWS_amp = [CSNSNR(isweep) 0];
    %         pfit.IWS_amp_jit = [0 CSNSNR(isweep)];
    %         pfit.SNRs = [0.2 CSNSNR(isweep) 0.4 -0.5];
    %     pfit.EPSP_amp = [0.25 1];
    %     pfit.CSN2AMN = CSNSNR(isweep);
    % [r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,BPfit,0,1);
    R = PB_scheme_setup();
    [R] = getFlaviesData(R);
    % pfit.EPSP_ampJit = [0 0 0];
    [xsims,dum,wflag] = PB_schema_simulate_REV_messaround(R,[],u,pfit,m);
    
    [dum,data_emp] = getFlaviesData(R);
    
    [~, feat_sim, wflag(2)] = R.obs.transFx(R,xsims);
    figure(1);  R.plot.outFeatFx(R,{R.data.feat_emp},{feat_sim},1);% PB_PlotOutputs
    
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
    [binStats_ampLag] = binDatabyRange(amp_sim,R.trans.amprange,lat_sim);
    binStats_ampLag(end,:) = [];
    
    plotBoundedScatter(binEdge2Mid(R.trans.amprange),binStats_ampLag(:,3)',binStats_ampLag(:,5)',amp_sim,lat_sim,cmap(1,:))
    
    xlabel('MEP Amplitude (z-scored)'); ylabel('Mean Corrected Latency')
    
    figure(4)
    %% Plot Phase Vs Amplitude
    subplot(3,2,1)
    [binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,amp_sim);
    binStats_phizAmp(end,:) = [];
    plotBoundedScatter(binEdge2Mid(R.trans.pirange),binStats_phizAmp(:,3)',binStats_phizAmp(:,5)',xsims.TMS_phase,amp_sim,cmap(1,:))
    xlim([-pi pi])
    title('Phase vs Amplitude'); xlabel('Pre-TMS Beta Phase'); ylabel('MEP Amplitude (z-scored)');
    %% Plot Phase Vs Latency
    subplot(3,2,3)
    hold on
    [binStats_phizLat] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,lat_sim);
    binStats_phizLat(end,:) = [];
    plotBoundedScatter(binEdge2Mid(R.trans.pirange),binStats_phizLat(:,3)',binStats_phizLat(:,5)',xsims.TMS_phase,lat_sim,cmap(3,:))
    xlim([-pi pi])
    title('Phase vs Latency'); xlabel('Pre-TMS Beta Phase'); ylabel('Mean Corrected Latency (ms)');
    %% Plot Phase Vs CoV
    subplot(3,2,5)
    s(1) = plot(feat_sim{2}(1,:),feat_sim{2}(2,:).*100); hold on
    %     s(2) = plot(R.data.feat_emp{2}(1,:),R.data.feat_emp{2}(2,:).*100,'k--');
    set(s(1),'LineWidth',2,'Color',cmap(1,:)); grid on
    
    s(3) = plot(feat_sim{4}(1,:),feat_sim{4}(2,:).*100);
    %     s(4) = plot(R.data.feat_emp{4}(1,:),R.data.feat_emp{4}(2,:).*100,'k--');
    set(s(3),'LineWidth',2,'Color',cmap(3,:)); grid on
    title('Phase vs Coefficient of Variation'); xlabel('Pre-TMS Beta Phase'); ylabel('% CoV')
    xlim([-pi pi])
    
    %% Plot Amplitude vs Latency
    subplot(3,2,2)
    [binStats_ampLag] = binDatabyRange(amp_sim,R.trans.amprange,lat_sim);
    binStats_ampLag(end,:) = [];
    plotBoundedScatter(binEdge2Mid(R.trans.amprange),binStats_ampLag(:,3)',binStats_ampLag(:,5)',amp_sim,lat_sim,cmap(6,:))
    xlim([-pi pi])
    title('Amplitude vs Latency'); xlabel('MEP Amplitude (z-scored)'); ylabel('Mean Corrected Latency')
    
    %% Plot Histograms
    subplot(3,2,4)
    h(1) = histogram(lat_sim(amp_sim>0),lat_supp,'Normalization','probability');
    h(1).FaceColor = cmap(7,:)
    hold on
    h(2) = histogram(lat_sim(amp_sim<0),lat_supp,'Normalization','probability');
    h(2).FaceColor = cmap(7,:).*[0.3 0.3 0.75];
    legend(h,{'Amp > 0','Amp < 0'})
    grid on
    xlabel('Mean Corrected Latency (ms)'); ylabel('Probability'); title('MEP Amplitude Distribution')
    xlim([-pi pi])
    set(gcf,'Position',[680    31   757   947]);
    
    %% Plot Model Outcomes
    figure(5)
    subplot(1,2,1)
    [binStats_phizCSNsync] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,xsims.CSN_sync);
    binStats_phizCSNsync(end,:) = [];
    plotBoundedScatter(binEdge2Mid(R.trans.pirange),binStats_phizCSNsync(:,3)',binStats_phizCSNsync(:,5)',xsims.TMS_phase,xsims.CSN_sync,cmap(4,:))
    title('Phase vs CSN Firing Synchrony'); xlabel('Pre-TMS Beta Phase'); ylabel('% CSN Sync')
    xlim([-pi pi])
    
    subplot(1,2,2)
    X = xsims.TMS_phase; Y = xsims.CSN_delay;
    Y = (Y./R.model.fsamp)*1000;
    X(Y<-30) = [];Y(Y<-30) = [];
    [binStats_phizCSN_delay] = binDatabyRange(X,R.trans.pirange,Y);
    binStats_phizCSN_delay(end,:) = [];
    plotBoundedScatter(binEdge2Mid(R.trans.pirange),binStats_phizCSN_delay(:,3)',binStats_phizCSN_delay(:,5)',X,Y,cmap(5,:))
    title('Phase vs CSN Firing Latency'); xlabel('Pre-TMS Beta Phase'); ylabel('CSN Spike Delay (ms)')
    xlim([-pi pi])
    set(gcf,'Position',[680   559   940   418]);
    
    figure(1)
end
