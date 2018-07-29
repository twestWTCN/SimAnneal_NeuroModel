%% Intracranial pressure - SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% Extrinsic Inputs arent signed - all summed together.
% Check extrinsic connections - something isnt being carried over from the
% DCM Prior structure
clear ; close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery'))
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Intracranial_pressure'))
% addpath('C:\spm12'); spm eeg; close all

%% Create 'data features' to be fitted
R = simannealsetup_ICP;
R.out.dag = '060718';
M = setupParVals_ICP; % values (model)
M.m = 1; % 1 "compartment"
P = setupParPriors_ICP_orig; % Prior means)
x = [M.Pic M.Ca];
M.Z       = {M.Pa M.dPadt M.vs M.I};
[xint T] = fx_simulateICP(R,x,[],P,M);
R.IntP.tvec = T';
R.data.feat_emp = xint;
R.data.feat_xscale = T;

clear xint T P

R.plot.save = 'True'
%% Try and infer the parameters of the generative model using simulated annealing
x = [M.Pic M.Ca];
P = setupParPriors_ICP;
[xobs1] = SimAn_ABC_110817(x,[],P,M,R);
gif_maker_siman(R)