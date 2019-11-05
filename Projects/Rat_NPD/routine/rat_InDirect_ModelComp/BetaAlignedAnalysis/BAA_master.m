%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
% simulateBurstData(R);
% load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim_F1.mat']); % This is the high density one
% load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat']); % This is the normal one
%% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R); % WILL ONLY WORK WITH FI
% %
% Lesion analysis
load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_MP.mat'],'MP')
BAA_sim_lesionExp(R,MP,32)

%% COMPUTE BETA BURSTS
R = computeBurstWrapper(R);
close all
%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
BB = AmpDurStatistics(R);

%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)

% This is a analysis of Beta inputs STR/M2
[R,MP] = BAA_sim_betaInputs(R,10,32);

% This is pertubation analysis to get PRCs- at the moment its not very
% clear what were looking for with it
[R] = BAA_sim_PRC(R,MP,500,1);
[R] = BAA_sim_fakeCloseLoop(R,MP,simtime,fresh);


% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
% 
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)




