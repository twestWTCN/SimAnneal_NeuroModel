function Qp = get_priors();
% setup priors
% I wave parameters
Qp.IWS(1,1) = 1; Qp.IWS(1,2) = 0.05; % Mean and std of Iwave Amp
Qp.IWS(2,1) = 1; Qp.IWS(2,2) = 0.05;

% EPSP decay parameters
Qp.EPSP_Tdecay = [0.005 0.005 0.005];
% EPSP sizes
Qp.EPSP_amp = [50 50 50];

% Spiking thresholds
Qp.SP_eps = [15 8];

