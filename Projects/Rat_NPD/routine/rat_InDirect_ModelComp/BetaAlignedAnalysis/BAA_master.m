%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
modID = 10;
%% Plot Model Fit
BAA_plotModelFit(R,modID,128);

% Lesion analysis
% BAA_sim_lesionExp_cortical_rhythms(R,modID,32,1)
BAA_sim_lesionExp(R,modID,32,1)


%% Plot Model Sweep Spectra
BAA_sim_ConnectionSweep_v2(R,modID,32,1)
% 'C:\Users\timot\Documents\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\InDrt_ModCompRev2\BetaBurstAnalysis\Data\BB_InDrt_ModCompRev2_ConnectionSweep_CON_1_feat_F1.mat'%
R = plotSweepSpectraWrapper(R); % You need the R produced here!
BAA_sim_ConnectionSweep_v2(R,modID,500,2) % Compute smaller list

%% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R); % WILL ONLY WORK WITH FI
R = plotSweepSpectraWrapper_M2_SI(R);
% %
simulateBurstData(R);

R = computeBurstWrapper(R);
%% Plot TimeLocked Statistics
BB.struccmap = linspecer(4);
TimeLockAnalysisMaster(R); %,BB,[15 17 20]); % [15 17 20] for STN_GPE [1 6 8]
OnsetEvolveAnalysisMaster(R)

% This is a analysis of Beta inputs STN/M2
[R,MP] = BAA_sim_betaInputs(R,10,32);

% This is pertubation analysis to get PRCs- at the moment its not very
% clear what were looking for with it
[R] = BAA_sim_PRC(R,MP,500,0);
[R] = BAA_sim_fakeCloseLoop(R,500,0);

[R] = BAA_sim_BetaPropagation(R,160,1); %remember reference of Boba and Hamacher 2015 for ZTE
 
 
% BB.range.RP = linspace(-pi,pi,7);
% BB = computeBetaBurstRPStats(R,BB);
% 
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)




