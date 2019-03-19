function [xsims,dum,wflag] = PB_schema_simulate2(R,~,~,p,m)
wflag = 0; dum = [];
fsamp = R.model.fsamp;

Qp = get_priors();
p.IWS = Qp.IWS .*exp(p.IWS);
p.EPSP_Tdecay = Qp.EPSP_Tdecay .*exp(p.EPSP_Tdecay);
p.EPSP_amp = Qp.EPSP_amp .*exp(p.EPSP_amp);
p.SP_eps = Qp.SP_eps .*exp(p.SP_eps);

%% Make Beta Signal (common to all CSNs)
% [t,tx_beta,px_beta] = makeBetaSignal(R.model.t_in);
% save('betaSignal','t','tx_beta','px_beta')
% load('betaSignal')
t = m.t;
tx_beta = m.tx_beta;
px_beta = m.px_beta;
%% Setup Kernels
% Make I waves
Iwave_win = linspace(0,0.05,(0.05*fsamp));
Iwave_amp = 15;
IwaveStat = p.IWS;
% TMS parameters
TMS_dur = floor(5*fsamp/1000);          % durations of TMS pulses (# time steps)
TMS_amp = 35;   % this matters only for plotting!
nPulses = R.model.nPulses;
[TMS_onsets,TMS_ind] = makeTMSSeries(t,nPulses,TMS_amp,TMS_dur,fsamp);
nTMS = numel(TMS_onsets);
% EPSPs
% For EPSP time constants of bi-exponential
% see: https://www.physiology.org/doi/pdf/10.1152/jn.00802.2005
epsp_win = linspace(0,0.075,(0.075*fsamp)); % EPSP time window
epsp_csn = makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(1)); % CSN PSP
epsp_amn = makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(2)); % AMN PSP
epsp_EMG = makeEPSP(epsp_win,0.002,p.EPSP_Tdecay(3)); % AMN PSP

%% Setup background activities
% Make I-wave "Sources"
IWS_n = 10;
tx_iws = [];
for cn = 1:IWS_n
    tx_iws(:,cn) = 1.5*randn(size(t));
    for i = 1:size(TMS_onsets,2)
        winInds = TMS_onsets(i):TMS_onsets(i)+length(Iwave_win)-1;
        if max(winInds)<size(tx_iws,1)
            Iwave = makeIwave(Iwave_win,Iwave_amp,IwaveStat);
            tx_iws(winInds,cn) = tx_iws(winInds,cn)' + Iwave;
        end
    end
end

% Make corticospinal neuron pool with noise + shared beta wave
CSN_n = 10;
tx_csn = [];
iws_thresh = p.SP_eps(1); % I wave thresh
for cn = 1:CSN_n
    [csn_P,nspike] = spike2psp(tx_iws(:,cn),iws_thresh,0.05*fsamp,p.EPSP_amp(1),epsp_csn);
    tx_csn(:,cn) = tx_beta + csn_P + randn(size(t)); % beta plus noise
    csn_scsrate(cn) = nspike/nTMS;
end

% Make alpha motor neuron pool with noise
AMN_n = 5;
srate = 100;
csn_thresh = p.SP_eps(2); % I wave thresh
tx_amn = [];
for an = 1:AMN_n
    % Random selection of layer above
    csn2amn = randperm(10,5);
    % spatial summation
    [amn_SumP,nspike] = spike2psp(tx_csn(:,csn2amn),csn_thresh,0.05*fsamp,p.EPSP_amp(2),epsp_amn);
    tx_amn(:,an) = amn_SumP + randn(size(t));
    amn_scsrate(an) = nspike/nTMS;
end


% Make EMG
amn_thresh = 4; % AMN threshold
while srate > 75
    amn_thresh = amn_thresh + 1;
    tx_EMG = [];
    [emg_SumP,nspike] = spike2psp(tx_amn,amn_thresh,0.05*fsamp,p.EPSP_amp(3),epsp_EMG);
    tx_EMG = emg_SumP + randn(size(t));
    tx_EMG = tx_EMG';
    EMG_scsrate = nspike/nTMS;
    srate = EMG_scsrate*100;
end

% Now loop through TMS pulse and measure MEP + delay
TMS_wind = ((100/1000)*fsamp); % window of action for TMS (to search for threshold)
for i = 1:size(TMS_onsets,2)-2
    windInds = TMS_onsets(i):TMS_onsets(i)+TMS_wind;
    TMS_phase(i) = px_beta(windInds(1));
    if max(windInds)<size(tx_EMG,1)
        wintx = tx_EMG(windInds);
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
MEP_srate = 100*(1-(sum(isnan(MEP_amp))./nTMS));
if MEP_srate<25 || MEP_srate>95
    wflag = 1;
    disp('Success Rate is too low!')
end
TMS_phase(MEP_onset<10) = [];
MEP_max(MEP_onset<10) = [];
MEP_amp(MEP_onset<10) = [];
MEP_onset(MEP_onset<10) = [];

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
xsims.tx_EMG = tx_EMG;

xsims.MEP_max = MEP_max;
xsims.MEP_mconset = MEP_mconset;
xsims.MEP_onset = MEP_onset;
xsims.MEP_amp = MEP_amp;

xsims.zMEP_max = zMEP_max;
xsims.zMEP_onset = zMEP_onset;
xsims.zMEP_amp = zMEP_amp;
plotTimeOutput(xsims,p)