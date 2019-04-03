function [R p m uc] = setup_model_testing(R)
%% MODEL 1 %%%
m.m = 1; % # of sources
m.x = [];
%% Prepare Priors
p.IWS_amp = [0 0];  % mean I1
p.IWS_amp_s = [0.15 0.15];

p.IWS_amp_jit = [0 0];  % mean I1
p.IWS_amp_jit_s = [0.15 0.15];


% EPSP decay parameters
p.EPSP_Tdecay = [0.20 0.32];
p.EPSP_Tdecay_s = [0.1 0.1];

% EPSP sizes
p.EPSP_amp = [0.6 0.33];
p.EPSP_amp_s = [1/4 1/4];

% EPSP Jitter
p.EPSP_ampJit = [0 0 ];
p.EPSP_ampJit_s = [1/4 1/4];

% Threshold for the success rates
p.SCRate = [0 0]; % CSN SR; AMN SR
p.SCRate_s = [1/8 1/8];

% SNR for beta; csn noise; amn noise; emg noise
p.SNRs = [0 0 0 0 ];
p.SNRs_s = [1/8 1/8 1/8 1/8];

p.CSN2AMN = 0;
p.CSN2AMN_s = 1/16;

p.CSN_n = 0;
p.CSN_n_s = 1/16;

p.AMN_n = 0;
p.AMN_n_s = 1/16;


% Spiking thresholds
% p.SP_eps = [0 0];
% p.SP_eps_s = [1/6 1/6];

p.obs.LF = 1;

uc= [];