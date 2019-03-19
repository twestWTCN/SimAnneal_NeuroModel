function m = setup_model()
m.m = 1;
% setup priors
% I wave parameters
m.IWS(1,1) = 1; Qp.IWS(1,2) = 0.05; % Mean and std of Iwave Amp
Qp.IWS(2,1) = 1; Qp.IWS(2,2) = 0.05;

% EPSP decay parameters
Qp.EPSP_Tdecay = [0.004 0.006 0.008];
% EPSP sizes
Qp.EPSP_amp = [100 10 25];

% Spiking thresholds
Qp.SP_eps = [15 8];

