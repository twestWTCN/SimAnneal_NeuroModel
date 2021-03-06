function [R p m uc] = setup_model_bpPriors(R)
%% MODEL 1 %%%
m.m = 1; % # of sources
m.x = [];
%% Prepare Priors
p.IWS(1,1) = 0.14;  % mean I1
p.IWS(1,2) = 0.45; % jitter I1 
p.IWS(2,1) = -0.48; % mean I2
p.IWS(2,2) = -0.31; % jitter I2

p.IWS_s(1,:) = [0.15 0.2];
p.IWS_s(2,:) = [0.15 0.17];

% EPSP decay parameters
p.EPSP_Tdecay = [0 0.20 0.32];
p.EPSP_Tdecay_s = [0.05 0.1 0.09];

% EPSP sizes
p.EPSP_amp = [1.21 0.6 0.33];
p.EPSP_amp_s = [0.13 0.46 0.22];

% EPSP Jitter
p.EPSP_ampJit = [0.17 0.56 -0.28];
p.EPSP_ampJit_s = [0.09 0.05 0.00];


% Spiking thresholds
% p.SP_eps = [0 0];
% p.SP_eps_s = [1/6 1/6];

p.obs.LF = 1;

uc= [];