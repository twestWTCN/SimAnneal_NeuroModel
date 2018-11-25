%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
% R.rootn = 'C:\Users\Tim\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.rootn = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'InDrt_ModComp';
daglist = {'NPD_InDrt_ModComp_M1','NPD_InDrt_ModComp_M2','NPD_InDrt_ModComp_M3',...
    'NPD_InDrt_ModComp_M4', 'NPD_InDrt_ModComp_M5', 'NPD_InDrt_ModComp_M6',...
    'NPD_InDrt_ModComp_M7', 'NPD_InDrt_ModComp_M8'};
mnum = 6; % best model
R.bandinits = {'\alpha','\beta_1','\beta_2'};
%% get data from model
 [R,permMod,xsimMod,permModHD,xsimModHD,permModSTR,xsimModSTR] = getSimData(R,daglist,mnum);
R.condname = {'Full','No HD','No M-STR'};
R.condcmap = linspecer(3);
R.cohband = 2;

save([R.rootn '\routine\rat_InDirect_ModelComp\BetaBurstAnalysis\Data\Indirect_Sims.mat'],...
    'R','permMod','xsimMod','permModHD','xsimModHD','permModSTR','xsimModSTR')
% OR
% load([R.rootn '\routine\rat_InDirect_ModelComp\BetaBurstAnalysis\Data\Indirect_Sims.mat'])
%% Get autospectra (for frequency finding)
AS = squeeze(permMod.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS_HD = squeeze(permModHD.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS_STR = squeeze(permModSTR.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
clear('permMod','permModHD','permModSTR')

BB = compute_BetaBursts_Simulated(R);


%% Do Beta Locked Analysis
PLVcmp = 0;
close all
F(1) = figure; F(2) = figure; 
BA_analysis(R,1,BB,xsimMod{1},AS_HD,PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model')

BA_analysis(R,2,BB,xsimModHD{1},AS_HD,PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

BA_analysis(R,3,BB,xsimModSTR{1},AS_STR,PLVcmp,F);
figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])

%% Investigate PLV High vs PLV Low
PLVcmp = 1;
close all
F(1) = figure; F(2) = figure;
BA_analysis(R,1,BB,xsimMod{1},AS_HD,PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model'); 

BA_analysis(R,2,BB,xsimModHD{1},AS_HD,PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

BA_analysis(R,3,BB,xsimModSTR{1},AS_STR,PLVcmp,F);
figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])

