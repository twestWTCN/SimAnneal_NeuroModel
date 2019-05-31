%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
% R = simulateBurstData(R);
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'])
%% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R)
% %
%% COMPUTE BETA BURSTS
[R,BB] = computeBurstWrapper(R);
close all
%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
BB = AmpDurStatistics(R,BB);

%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R,BB,[5 4 6])
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
%
[R,MP] = BAA_sim_betaInputs(R,10,32);
[R] = BAA_sim_PRC(R,MP,500);