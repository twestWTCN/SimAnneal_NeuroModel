function [xsims,dum,wflag] = PB_schema_simulate_REV(R,~,uc,p,m)
wflag = 0; dum = [];
fsamp = R.model.fsamp;

Qp = get_priors();
p.IWS = Qp.IWS .*exp(p.IWS);
p.EPSP_Tdecay = Qp.EPSP_Tdecay .*exp(p.EPSP_Tdecay);
p.EPSP_amp = Qp.EPSP_amp .*exp(p.EPSP_amp);
p.EPSP_ampJit = Qp.EPSP_ampJit .*exp(p.EPSP_ampJit);

p.SCRate = Qp.SCRate .*exp(p.SCRate);

% p.SP_eps = Qp.SP_eps .*exp(p.SP_eps);

%% Make Beta Signal (common to all CSNs)
% [t,tx_beta,px_beta] = makeBetaSignal(R.model.t_in);
% save('betaSignal','t','tx_beta','px_beta')
% load('betaSignal')
t = m.t;
tx_beta = 5.*m.tx_beta;
px_beta = m.px_beta;
%% Setup Kernels
% Make I waves
Iwave_win = linspace(0,0.05,(0.05*fsamp));
Iwave_amp = 25;
IwaveStat = p.IWS;
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
epsp_csn = p.EPSP_amp(1).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(1)); % CSN PSP
epsp_amn = p.EPSP_amp(2).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(2)); % AMN PSP
epsp_EMG = p.EPSP_amp(3).*makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(3)); % AMN PSP

%% LEVEL 1: Make corticospinal neuron pool with noise + shared beta wave
% Setup the I-wave depolarizations
CSN_n = 10;
tx_csn = nan(size(t,2),CSN_n);
for cn = 1:CSN_n
    tx_csn(:,cn) = tx_beta + 1.5*randn(size(t));
    for i = 1:size(TMS_onsets,2)
    winInds = TMS_onsets(i):TMS_onsets(i)+length(Iwave_win)-1;
        if max(winInds)<size(tx_csn,1)
            Iwave = makeIwave(Iwave_win,Iwave_amp,IwaveStat);
            tx_csn(winInds,cn) = tx_csn(winInds,cn)' + Iwave;
        end
    end
end

srate = 0;
csn_thresh = 0.5.*max(tx_csn(:)); % I wave thresh
i = 0;
% Get spike firings of the CSN pool responding to CSN Input
while (srate < 30  || srate > p.SCRate(1)) && csn_thresh < max(tx_csn(:))
    i = i + 1;
    csn_thresh = csn_thresh + (csn_thresh./50); % find upper-bound!
    in_spT = []; csn_scsrate = [];
    for cn = 1:CSN_n
        [csn_spT{cn},nspike,nspike_inWin] = findSpike(tx_csn(:,cn),csn_thresh,0.075*fsamp,TMS_win);
        csn_scsrate(cn) = 100.*(nspike_inWin/nTMS);
    end
    srate = mean(csn_scsrate);
end


%% LEVEL 2: Make alpha motor neuron pool with noise
AMN_n = 5;
% Now convolve the spike times with the EPSP kernel to get the
% depolarizations of the AMN
tx_amn = nan(size(t,2),AMN_n);
for cn = 1:AMN_n
    csn2amn = randperm(CSN_n,4);
    sum_spT = sort(vertcat(csn_spT{csn2amn})); % summated spike times for AMN 
    csn_P = convSpikePSP(sum_spT,epsp_amn,p.EPSP_ampJit(1),t); %tx_iws(:,cn));
    tx_amn(:,cn) = (csn_P + randn(size(t)))'; % beta plus noise
end

srate = 0;
amn_thresh = 0.5.*max(tx_amn(:)); % I wave thresh
i = 0;
% Get the AMN Spike Times
while (srate < 15  || srate > 65) && amn_thresh < max(tx_amn(:))
    i = i + 1;
    amn_thresh = amn_thresh + (csn_thresh./50);
    in_spT = []; amn_scsrate = [];
    for an = 1:AMN_n
        % spatial summation
        [amn_spT{an},nspike,nspike_inWin] = findSpike(tx_amn(:,an),amn_thresh,0.075*fsamp,TMS_win);
        amn_scsrate(an) = 100.*(nspike_inWin/nTMS);
    end
    srate = max(amn_scsrate);
end

%% Level 3: Get the EMG depolarization (MEP!)
% Now convolve the spike times with the EPSP kernel
tx_emg = nan(size(t,2),1);
amn2emg = 1:AMN_n;
sum_spT = sort(vertcat(amn_spT{amn2emg})); % summated spike times for AMN
emg_P = convSpikePSP(sum_spT,epsp_EMG,p.EPSP_ampJit(1),t); %tx_iws(:,cn));
tx_emg = (emg_P + randn(size(t)))'; % beta plus noise


%% MEP ANALYSIS
% Now loop through TMS pulse and measure MEP + delay
TMS_wind = ((100/1000)*fsamp); % window of action for TMS (to search for threshold)
for i = 1:size(TMS_onsets,2)-2
    windInds = TMS_onsets(i):TMS_onsets(i)+TMS_wind;
    TMS_phase(i) = px_beta(windInds(1));
    if max(windInds)<size(tx_emg,1)
        wintx = tx_emg(windInds);
        wintx(wintx<5) = 0;
    else
        wintx = 0;
    end
    if any(wintx)
        MEP_max(i) =  1000*(find(wintx == max(wintx),1,'first')./fsamp);
        MEP_onset(i) =  1000*(find(wintx > 0,1,'first')./fsamp);
        MEP_amp(i) = max(wintx);
    else % failed to ellict MEP
        MEP_max(i) = NaN;
        MEP_onset(i) =  NaN;
        MEP_amp(i) = NaN;
    end
end
tx_emg(tx_emg<5) = 0;
MEP_srate = 100*(sum(~isnan(MEP_amp))./nTMS);
if MEP_srate<15 || MEP_srate>85
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



plot(t,tx_csn,'r'); hold on; plot(t,tx_amn,'b'); hold on; plot(t,tx_emg,'k')


plotTimeOutput(xsims,p)