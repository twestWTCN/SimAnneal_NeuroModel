function [R p m uc] = setup_model(R)
%% MODEL 1 %%%
m.m = 1; % # of sources
m.x = [];
%% Prepare Priors
p.IWS(1,1) = 0;  % mean I1
p.IWS(1,2) = 0; % jitter I1 
p.IWS(2,1) = 0; % mean I2
p.IWS(2,2) = 0; % jitter I2

p.IWS_s(1,:) = [1/8 1/8];
p.IWS_s(2,:) = [1/8 1/8];

% EPSP decay parameters
p.EPSP_Tdecay = [0 0 0];
p.EPSP_Tdecay_s = [1/8 1/8 1/8];

% EPSP sizes
p.EPSP_amp = [0 0 0];
p.EPSP_amp_s = [1/4 1/4 1/4];

% Spiking thresholds
p.SP_eps = [0 0];
p.SP_eps_s = [1/6 1/6];

p.obs.LF = 1;

uc= [];