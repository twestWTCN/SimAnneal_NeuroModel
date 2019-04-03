function Qp = get_priors()
% setup priors
% I wave parameters
Qp.IWS_amp = [1 0.75];  % Mean of 1st and 2nd I wave
Qp.IWS_amp_jit = [0.05 0.05];  % Mean of 1st and 2nd I wave

% EPSP decay parameters
Qp.EPSP_Tdecay = [0.008 0.008]; % AMN; CSN
% EPSP sizes
Qp.EPSP_amp = [50 50]; % AMN; CSN
% EPSP jitter
Qp.EPSP_ampJit = [0.5 0.5]; % CSN; AMN
% Spiking thresholds
Qp.SCRate = [65 60]; % success rate
% SNR
Qp.SNRs = [5 2 2 2]; % beta; CSN; AMN; EMG

Qp.CSN_n = 10;
Qp.AMN_n = 8;
Qp.CSN2AMN = 5;
% MASON 1991 - EPSP variability 