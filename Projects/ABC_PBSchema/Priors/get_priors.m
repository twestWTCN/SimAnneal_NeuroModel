function Qp = get_priors()
% setup priors
% I wave parameters
Qp.IWS_amp = [1 0.8];  % Mean of 1st and 2nd I wave
Qp.IWS_amp_jit = [0.1 0.1];  % Mean of 1st and 2nd I wave

% EPSP decay parameters
Qp.EPSP_Tdecay = [0.008 0.003]; % AMN; CSN
% EPSP sizes
Qp.EPSP_amp = [35 35]; % AMN; CSN
% EPSP jitter
Qp.EPSP_ampJit = [0.1 0.1]; % CSN; AMN
% Spiking thresholds
Qp.SCRate = [85 75]; % success rate
% SNR
Qp.SNRs = [2 2 2 2]; % beta; CSN; AMN; EMG

Qp.CSN_n = 10;
Qp.AMN_n = 8;
Qp.CSN2AMN = 4;
% MASON 1991 - EPSP variability 