function Qp = get_priors();
% setup priors
% I wave parameters
Qp.IWS(1,1) = 1; Qp.IWS(1,2) = 0.05; % Mean and std of Iwave Amp
Qp.IWS(2,1) = 1; Qp.IWS(2,2) = 0.05;

% EPSP decay parameters
Qp.EPSP_Tdecay = [0.015 0.008 0.008];
% EPSP sizes
Qp.EPSP_amp = [50 50 50];
% EPSP jitter
Qp.EPSP_ampJit = [0.1 0.1 0.1];
% Spiking thresholds
Qp.SP_eps = [15 8];

Qp.SCRate = [65 60 75];


% MASON 1991 - EPSP variability 