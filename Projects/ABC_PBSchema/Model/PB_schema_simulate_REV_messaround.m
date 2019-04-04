function [xsims,dum,wflag] = PB_schema_simulate_REV_messaround(R,~,uc,p,m)
wflag = 0; dum = [];
fsamp = R.model.fsamp;
rng(2121)
Qp = get_priors();
p.IWS_amp = Qp.IWS_amp .*exp(p.IWS_amp);
p.IWS_amp_jit = Qp.IWS_amp_jit .*exp(p.IWS_amp_jit);
p.EPSP_Tdecay = Qp.EPSP_Tdecay .*exp(p.EPSP_Tdecay); % EPSP decay rate
p.EPSP_amp = Qp.EPSP_amp .*exp(p.EPSP_amp); % EPSP amplitudes AMN; MEP
p.EPSP_ampJit = Qp.EPSP_ampJit .*exp(p.EPSP_ampJit); % EPSP variability for csn; amn
p.SNRs = Qp.SNRs.*exp(p.SNRs); % Scaling for background: beta, CSN noise, AMN noise, EMG noise
p.SCRate = Qp.SCRate .*exp(p.SCRate); % Success rates thresholds (CSNs, AMNs)
p.AMN_n = round(Qp.AMN_n .*exp(p.AMN_n));
p.CSN_n = round(Qp.CSN_n .*exp(p.CSN_n));
p.CSN2AMN = round(Qp.CSN2AMN .*exp(p.CSN2AMN));

if p.AMN_n>10
    p.AMN_n = 10;
end
if p.AMN_n <= 0
    p.AMN_n = 1;
end
if p.CSN_n>10
    p.CSN_n = 10;
end
if p.CSN_n <= 0
    p.CSN_n = 1;
end
% p.SP_eps = Qp.SP_eps .*exp(p.SP_eps);
% MY CHOICES
p.SCRate  = [85 75];
p.IWS_amp = [1 0.8];
p.IWS_amp_jit = [0.1 0.1];
p.CSN_n = 10;
p.CSN2AMN = 4;
p.EPSP_Tdecay(2) = 0.003;
p.EPSP_ampJit = [0.1 0.1];
p.SNRs(1) = 2;
% p.IWS_amp = [2.7 1.5];
load('C:\Users\timot\Documents\GitHub\SimAnneal_NeuroModel\Projects\ABC_PBSchema\testing\sample_save\working_PBV2.mat')

%% Make Beta Signal (common to all CSNs)
% [t,tx_beta,px_beta] = makeBetaSignal(R.model.t_in);
% save('betaSignal','t','tx_beta','px_beta')
% load('betaSignal')
t = m.t;

tx_beta = p.SNRs(1).*m.tx_beta;
px_beta = m.px_beta;
%% Setup Kernels
% Make I waves
Iwave_win = linspace(0,0.05,(0.05*fsamp));
Iwave_amp = 25; % absolute I wave amp
IwaveStat = [p.IWS_amp;p.IWS_amp_jit]';
% TMS parameters
% TMS_dur = floor(5*fsamp/1000);          % durations of TMS pulses (# time steps)
% TMS_winL= floor(100*fsamp/1000);          % durations of TMS pulses (# time steps)
%
% TMS_amp = 35;   % this matters only for plotting!
% nPulses = R.model.nPulses;
% [TMS_onsets,TMS_ind,TMS_win] = makeTMSSeries(t,nPulses,TMS_amp,TMS_dur,TMS_winL,fsamp);
% nTMS = numel(TMS_onsets);
TMS_onsets = uc.TMS_onsets;
TMS_ind = uc.TMS_ind;
TMS_win = uc.TMS_win;
nTMS = uc.nTMS;
% EPSPs
% For EPSP time constants of bi-exponential
% see: https://www.physiology.org/doi/pdf/10.1152/jn.00802.2005
epsp_win = linspace(0,0.075,(0.075*fsamp)); % EPSP time window
% epsp_csn = p.EPSP_amp(1).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(1)); % CSN PSP
epsp_amn = p.EPSP_amp(1).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(1)); % AMN PSP
epsp_EMG = p.EPSP_amp(2).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(2)); % EMG PSP

%% LEVEL 1: Make corticospinal neuron pool with noise + shared beta wave
% Setup the I-wave depolarizations
tx_csn = nan(size(t,2),p.CSN_n);
for cn = 1:p.CSN_n
    tx_csn(:,cn) = tx_beta + p.SNRs(2)*randn(size(t));
    for i = 1:size(TMS_onsets,2)
        winInds = TMS_onsets(i):TMS_onsets(i)+length(Iwave_win)-1;
        if max(winInds)<size(tx_csn,1)
            Iwave = makeIwave(Iwave_win,Iwave_amp,IwaveStat);
            tx_csn(winInds,cn) = tx_csn(winInds,cn)' + Iwave;
        end
    end
end
a = plot(t,tx_csn,'r');
hold on
xleg{1} = 'CSN Potential';
pst(1) = a(1);
srate = 0;
csn_thresh = 1.05.*max(tx_beta); % I wave thresh
i = 0;
% Get spike firings of the CSN pool responding to CSN Input
while (srate < 30  || srate > p.SCRate(1)) && csn_thresh < max(tx_csn(:))
    i = i + 1;
    csn_thresh = csn_thresh + (csn_thresh./50); % find upper-bound!
    in_spT = []; csn_scsrate = [];
    for cn = 1:p.CSN_n
        [csn_spT{cn},nspike,nspike_inWin] = findSpike(tx_csn(:,cn),csn_thresh,0.075*fsamp,TMS_win,t);
        csn_scsrate(cn) = 100.*(nspike_inWin/nTMS);
    end
    srate = mean(csn_scsrate);
end

%% LEVEL 2: Make alpha motor neuron pool with noise
% Now convolve the spike times with the EPSP kernel to get the
% depolarizations of the AMN
tx_amn = nan(size(t,2),p.AMN_n); spT = [];
for cn = 1:p.AMN_n
    csn2amn = randperm(p.CSN_n,p.CSN2AMN );
    sum_spT = sort(vertcat(csn_spT{csn2amn})); % summated spike times for AMN
    [csn_P,spT(:,cn)] = convSpikePSP(sum_spT,epsp_amn,p.EPSP_ampJit(1),t); %tx_iws(:,cn));
    tx_amn(:,cn) = (csn_P + p.SNRs(3).*randn(size(t)))'; % beta plus noise
end
hold on
spT(spT>0) = 1;
a = plot(t,spT.*csn_thresh,'r');
pst(2)  = a(1);
a = plot(t,tx_amn(:,:),'b');
pst(3)  = a(1);
xleg{2} = 'CSN Firing';
xleg{3} = 'AMN Potential';

srate = 0;
amn_thresh = 0.25.*max(tx_amn(:)); % I wave thresh
i = 0;
% Get the AMN Spike Times
while (srate < 15  || srate > p.SCRate(2)) && amn_thresh < max(tx_amn(:))
    i = i + 1;
    amn_thresh = amn_thresh + (amn_thresh./50);
    in_spT = []; amn_scsrate = [];
    for an = 1:p.AMN_n
        % spatial summation
        [amn_spT{an},nspike,nspike_inWin] = findSpike(tx_amn(:,an),amn_thresh,0.075*fsamp,TMS_win);
        amn_scsrate(an) = 100.*(nspike_inWin/nTMS);
    end
    srate = max(amn_scsrate);
end

%% Level 3: Get the EMG depolarization (MEP!)
% Now convolve the spike times with the EPSP kernel
tx_emg = nan(size(t,2),1);
amn2emg = 1:p.AMN_n;
sum_spT = sort(vertcat(amn_spT{amn2emg})); % summated spike times for AMN
[emg_P,spT] = convSpikePSP(sum_spT,epsp_EMG,p.EPSP_ampJit(2),t); %tx_iws(:,cn));
tx_emg = (emg_P +  p.SNRs(4).*randn(size(t)))'; % beta plus noise

spT(spT>0) = 1;
a = plot(t,spT.*amn_thresh,'b');
pst(4) = a(1);
xleg{4} = 'AMN Firing';

pst(5) = plot(t,tx_emg,'g');
xleg{5} = 'MEP';
legend(pst,xleg)
%% MEP ANALYSIS
% Now loop through TMS pulse and measure MEP + delay
TMS_wind = ((100/1000)*fsamp); % window of action for TMS (to search for threshold)
MEP_srate = 0;
emg_thresh = 0.9*max(epsp_EMG); %0.1.*max(tx_emg(:)); % I wave thresh
L = 0;
% while (MEP_srate < 45  || MEP_srate > p.SCRate(3)) && emg_thresh < max(tx_emg(:))
%     L = L + 1;
emg_thresh = emg_thresh; % + (emg_thresh./50);
for i = 1:size(TMS_onsets,2)-2
    windInds = TMS_onsets(i):TMS_onsets(i)+TMS_wind;
    TMS_phase(i) = px_beta(windInds(1));
    [spSync, spFire, spDelay] = computeLayerSync(windInds,csn_spT);
    if max(windInds)<size(tx_emg,1)
        wintx = tx_emg(windInds);
        wintx(wintx<emg_thresh) = 0;
    else
        wintx = 0;
    end
    if any(wintx)
        CSN_sync(i) = spSync*100;
        CSN_delay(i) = nanmean(spDelay);
        
        MEP_max(i) =  1000*(find(wintx == max(wintx),1,'first')./fsamp);
        MEP_onset(i) =  1000*(find(wintx > 0,1,'first')./fsamp);
        MEP_amp(i) = max(wintx);
    else % failed to ellict MEP
        CSN_sync(i) = NaN;
        CSN_delay(i) = NaN;
        
        MEP_max(i) = NaN;
        MEP_onset(i) =  NaN;
        MEP_amp(i) = NaN;
    end
end
MEP_srate = 100*(sum(~isnan(MEP_amp))./nTMS);
% end
tx_emg(tx_emg<emg_thresh) = 0;

if MEP_srate<45 || MEP_srate>85
    wflag = 1;
    disp('Success Rate is too low!')
end
TMS_phase(MEP_onset<10) = [];
MEP_max(MEP_onset<10) = [];
MEP_amp(MEP_onset<10) = [];
MEP_onset(MEP_onset<10) = [];


TMS_phase(MEP_onset>60) = [];
MEP_max(MEP_onset>60) = [];
MEP_amp(MEP_onset>60) = [];
MEP_onset(MEP_onset>60) = [];


MEP_max = MEP_max;
MEP_onset = MEP_onset;
MEP_mconset = (MEP_onset-nanmean(MEP_onset));
MEP_amp = MEP_amp;

zMEP_max = (MEP_max-nanmean(MEP_max))./nanstd(MEP_max);
zMEP_onset = (MEP_onset-nanmean(MEP_onset))./nanstd(MEP_onset);
zMEP_amp = (MEP_amp-nanmean(MEP_amp))./nanstd(MEP_amp);

xsims.t = t;
xsims.TMS_phase = TMS_phase;
xsims.TMS_ind = TMS_ind;
xsims.tx_csn = tx_csn;
xsims.tx_amn = tx_amn;
xsims.tx_emg = tx_emg;

xsims.MEP_max = MEP_max;
xsims.MEP_mconset = MEP_mconset;
xsims.MEP_onset = MEP_onset;
xsims.MEP_amp = MEP_amp;

xsims.zMEP_max = zMEP_max;
xsims.zMEP_onset = zMEP_onset;
xsims.zMEP_amp = zMEP_amp;
xsims.csn_thresh = csn_thresh;
xsims.amn_thresh = amn_thresh;
xsims.CSN_sync = CSN_sync;
xsims.CSN_delay = CSN_delay;

% scatter(TMS_phase,MEP_amp); yyaxis right; scatter(TMS_phase,MEP_mconset)

% xedg = -3:0.3:3;
% histogram(zMEP_onset(zMEP_amp>1.2),xedg,'Normalization','probability');
% hold on; histogram(zMEP_onset(zMEP_amp<1.2),xedg,'Normalization','probability')
% plot(t,tx_csn,'r'); hold on; plot(t,tx_amn,'b'); hold on; plot(t,tx_emg,'k')
% plotTimeOutput(xsims,p)