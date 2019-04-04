function precomputeTMS(R,m)
t = makeBetaSignal(R.model.t_in);
fsamp = R.model.fsamp;

% TMS parameters
TMS_dur = floor((1*fsamp)/1000);          % durations of TMS pulses (# time steps)
TMS_winL= floor((120*fsamp)/1000);          % durations of TMS pulse window (# time steps)

TMS_amp = 35;   % this matters only for plotting!
nPulses = R.model.nPulses;

for i = 1:25
    [TMS_onsets,TMS_ind,TMS_win] = makeTMSSeries(t,nPulses,TMS_amp,TMS_dur,TMS_winL,fsamp);
    nTMS = numel(TMS_onsets);
    u = [];
    u.TMS_onsets = TMS_onsets;
    u.TMS_ind = TMS_ind;
    u.TMS_win = TMS_win;
    u.nTMS = nTMS;
    saveMkPath([R.rootn '\Inputs\TMSsig\TMSsignal_' num2str(i)],u)
end