%% Master - Compute Beta Burst Aligned Analysis for Indirect Comp %%%
clear; close all
% simAnnealAddPaths()

%% Setup
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condname = {'Full','No HD','No M-STR'};
R.condcmap = linspecer(3);
R.cohband = 2;

modlist = [6 12 15];
i = 0;
for modID = modlist
    i = i + 1;
    R.out.tag = 'InDrt_ModComp';
    R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
    [R,m,p,parBank] = loadABCData(R);
    %% get data from model
    R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
    [R,permMod(i),xsimMod(i)] = getSimData_v2(R,modID,500);
end
R.out.tagOld = 'rat_InDirect_ModelComp';
mkdir([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data'])
save([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],...
    'permMod','xsimMod')
% OR
% load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'])
%% Get autospectra (for frequency finding)
AS{1} = squeeze(permMod{1}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS{2} = squeeze(permMod{2}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS{3} = squeeze(permMod{3}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
clear('permMod','permModHD','permModSTR')

%% CROSS OVER
R.BBA_path = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\BetaBursts';
addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\statistics');
addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Plotting');
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\superbar'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\PEB_DFA');
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MLDFA');

R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condname = {'Full','No HD','No M-STR'};
R.condcmap = linspecer(3);
R.cohband = 2;
%%
BB = compute_BetaBursts_Simulated(R);

%% Do Beta Locked Analysis
addpath(R.BBA_path)

PLVcmp = 0;
close all
F(1) = figure; F(2) = figure;
BA_analysis(R,1,BB,xsimMod{1}{1},AS{1},PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model')

BA_analysis(R,2,BB,xsimMod{2}{1},AS{2},PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

BA_analysis(R,3,BB,xsimMod{3}{1},AS{3},PLVcmp,F);
figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])

%% Investigate PLV High vs PLV Low
PLVcmp = 1;
close all
F(1) = figure; F(2) = figure;
BA_analysis(R,1,BB,xsimMod{1}{1},AS{1},PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model');

BA_analysis(R,2,BB,xsimMod{2}{1},AS{2},PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

BA_analysis(R,3,BB,xsimMod{3}{1},AS{3},PLVcmp,F);
figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])

