%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
R = simulateBurstData(R);
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim_F1.mat']); % This is the high density one
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat']); % This is the normal one
%% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R); % WILL ONLY WORK WITH FI
% %
%% COMPUTE BETA BURSTS
[R,BB] = computeBurstWrapper(R);
close all
%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
BB = AmpDurStatistics(R,BB);

BB.range.RP = linspace(-pi,pi,7);
BB = computeBetaBurstRPStats(R,BB);
%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R,BB,[14 16 20])
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
%
[R,MP] = BAA_sim_betaInputs(R,10,32);
[R] = BAA_sim_PRC(R,MP,500);


BAA_sim_lesionExp(R,MP,32)