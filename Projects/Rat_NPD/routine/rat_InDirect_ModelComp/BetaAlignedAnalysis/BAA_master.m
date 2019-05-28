%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Simulate Data
R = simulateBurstData(R);
% % load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep.mat'])
% % 
% %% Plot Model Sweep Spectra
plotSweepSpectraWrapper(R)
% %  
%% COMPUTE BETA BURSTS
[R,BB] = computeBurstWrapper(R);

%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
figureBB = AmpDurStatistics(R,BB);

%% Plot TimeLocked Statistics
% BB.struccmap = linspecer(4);
% % Get Autospectra
% for i = 1:4;        AS{1}(i,:) = squeeze(feat{1}{1}(1,i,i,1,:));    end
% for cond = 2:numel(R.condname)
%     for i = 1:4;        AS{cond}(i,:) = squeeze(feat{cond}(1,i,i,1,:));    end
% end
% 
% SimulationTLBurstAnaly(R,BB,[4 1 5],AS)
% % TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)
% 
