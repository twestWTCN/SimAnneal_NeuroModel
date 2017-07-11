%% Intracranial pressure - SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% Extrinsic Inputs arent signed - all summed together.
% Check extrinsic connections - something isnt being carried over from the
% DCM Prior structure
clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))


%% Create 'data features' to be fitted
R = simannealsetup_ICP;
M = setupParVals_ICP; % values (model)
M.m = 1; % 1 "compartment"
P = setupParPriors_ICP_orig; % Prior means)
x = [M.Pic M.Ca];
M.Z       = {M.Pa M.dPadt M.vs M.I};
[xint T] = fx_simulateICP(x,M,P);
R.data.feat_emp = xint;
R.data.feat_xscale = T;

clear xint T P


%% Try and infer the parameters of the generative model using simulated annealing
x = [M.Pic M.Ca];
P = setupParPriors_ICP;
[xobs1] = SimAn_100717(x,[],P,M,R);
% gif_maker_siman(R)