%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
R = simulateBurstData(R);
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep.mat'])
% 
%% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R)
% %  
%% COMPUTE BETA BURSTS
[R,BB] = computeBurstWrapper(R);
close all
%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
BB = AmpDurStatistics(R,BB);

%% Plot TimeLocked Statistics
    BB.struccmap = linspecer(4);
    TimeLockAnalysisMaster(R,BB,[2 1 3],[])
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
% 
 [R] = BAA_sim_betaInputs(R,10,25)