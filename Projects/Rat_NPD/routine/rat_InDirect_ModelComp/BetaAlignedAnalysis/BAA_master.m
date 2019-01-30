%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()
% addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\statistics');
% addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Plotting');

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condcmap = linspecer(3);
R.cohband = 2;

%% Simulate Data
R = BAA_sim_discreteModels(R,12,250);
% OR load
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],'xsimMod','feat')    
 
% Setup condition variables
R.condname = {'Full','1% M2->STN','150% M2->STN','1% STR->GPe','150% STR->GPe','Empirical'};
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 6 8 16 18],:);
R.condcmap(6,:) = [0 0 0];
R = BAA_sim_ConnectionSweep(R,12,32);

%% COMPUTE BETA BURSTS
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
R.condname = {'Full','1% M2->STN','150% M2->STN','1% STR->GPe','150% STR->GPe','Empirical'};
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 6 8 16 18],:);
R.condcmap(6,:) = [0 0 0];

[R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
BB = compute_BurstThreshold(R,BB,0);
% Compute Time Decomposition
% Compute Dur/AMP Statistics
close all

%% Plot Model Sweep Spectra
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep.mat'],'feat_HD','feat_STR_GPe')
figure
subplot(1,2,1)
 plotSweepSpectra(R.frqz,feat_HD,feat,cmap([4 7],:),R.condcmap(1,:),{R.condname{[1 2 3 6]}})
hold on
subplot(1,2,2)
 plotSweepSpectra(R.frqz,feat_STR_GPe,feat,cmap([15 18],:),R.condcmap(1,:),{R.condname{[1 4 5 6]}})
set(gcf,'Position',[684         643        1024         224])

%% Plot Burst Statistics
load([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
figure
BB = AmpDurStatistics(R,BB);

%% Plot TimeLocked Statistics
 BB.struccmap = linspecer(4);
% Get Autospectra
for i = 1:4;        AS{1}(i,:) = squeeze(feat{1}{1}(1,i,i,1,:));    end
for cond = 2:numel(R.condname)
for i = 1:4;        AS{cond}(i,:) = squeeze(feat{cond}(1,i,i,1,:));    end
end

SimulationTLBurstAnaly(R,BB,[4 1 5],AS)
% TimeLockedBetaPLVAnalysis(R,BB,xsimMod,AS)

