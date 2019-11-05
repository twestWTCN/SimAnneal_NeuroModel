%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
simulateBurstData(R);

% Lesion analysis
BAA_sim_lesionExp(R,32)

%% Plot Model Sweep Spectra
R = plotSweepSpectraWrapper(R); % You need the R produced here!
BAA_sim_ConnectionSweep_v2(R,10,500,2) % Compute smaller list

R = computeBurstWrapper(R);
%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)

% This is a analysis of Beta inputs STN/M2
[R,MP] = BAA_sim_betaInputs(R,10,32);

% This is pertubation analysis to get PRCs- at the moment its not very
% clear what were looking for with it
[R] = BAA_sim_PRC(R,MP,500,1);
[R] = BAA_sim_fakeCloseLoop(R,500,0);

[R] = BAA_sim_BetaPropagation(R,120,1); %remember reference of Boba and Hamacher 2015 for ZTE
 
 
% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
% 
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)




