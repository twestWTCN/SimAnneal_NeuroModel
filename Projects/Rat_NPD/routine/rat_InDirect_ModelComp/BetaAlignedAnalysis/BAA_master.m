%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()
% addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\statistics');
% addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Plotting');

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condname = {'Full','No HD','No STR-GPe'};
R.condcmap = linspecer(3);
R.cohband = 2;

%% Simulate Data
R = BAA_sim_data(R,12,128);
% OR load
R.out.tagOld = 'rat_InDirect_ModelComp'; % This is hack for old naming system
load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],...
    'permMod','xsimMod')
 
% Setup condition variables
R.condname = {'Full (M12)','No HD (M6)','No M-STR (M18)'};
NPDsel = [12 6 18];
cmap = brewermap(18,'Spectral');
R.condcmap = cmap(NPDsel,:);

%% COMPUTE BETA BURSTS
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
[R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
% Compute Time Decomposition
% Compute Dur/AMP Statistics
% Get Autospectra
for i = 1:4;        AS{1}(i,:) = squeeze(permMod{1}.feat_rep{1}(1,i,i,1,:));    end
for i = 1:4;        AS{2}(i,:) = squeeze(permMod{2}.feat_rep{1}(1,i,i,1,:));    end
for i = 1:4;        AS{3}(i,:) = squeeze(permMod{3}.feat_rep{1}(1,i,i,1,:));    end

load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
close all
BB = AmpDurStatistics(R,BB);

BB.struccmap = linspecer(4);

SimulationTLBurstAnaly(R,BB,xsimMod,AS)
% TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)

