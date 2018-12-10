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

i = 0;
i = i + 1;
modID = 12;
R.out.tag = 'InDrt_ModComp';
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
%% get data from model
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,permMod,xsimMod] = getSimData_v2(R,8,20);

close all
clear r2mean feat_sim xsims dfx
Alist = 0.15:0.15:0.9; % [0,0.25,0.33,0.5,1];
Alist(3) = 0.4483;
cmap = linspecer(numel(Alist));
for i = 1:length(Alist)
    pnew = permMod{1}.par_rep{1};
    pnew.A{1}(4,1) = Alist(i);
    [r2mean{i},dum,feat_sim{i},xsims{i}] = computeSimData(R,m,pnew,2000,0); 
    fx(1,:) = squeeze(feat_sim{i}(1,1,1,1,:));
    fx(2,:) = squeeze(feat_sim{i}(1,4,4,1,:));
    fx(3:4,:) = squeeze(feat_sim{i}(1,1,4,2:3,:));
    dfx{i} = fx;
    bp(i) = sum(fx(1,R.frqz>14 & R.frqz<24));
    for p = 1:4
        subplot(1,4,p)
        plot(R.frqz,fx(p,:),'color',cmap(i,:));
        hold on
    end
end  

    subplot(1,4,1); ylim([0 18]); xlabel('Frequency Hz'); ylabel('Amplitude'); title('M2 Autospectra')
    subplot(1,4,2); ylim([0 18]); xlabel('Frequency Hz'); ylabel('Amplitude'); title('STN Autospectra')
    subplot(1,4,3); ylim([0 0.2]); xlabel('Frequency Hz'); ylabel('NPD'); title('M2 -> STN NPD')
    subplot(1,4,4); ylim([0 0.2]); xlabel('Frequency Hz'); ylabel('NPD'); title('STN -> M2 NPD')
    
R.out.tagOld = 'rat_InDirect_ModelComp';
mkdir([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data'])
save([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],...
    'permMod','xsimMod')
% OR
% load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'])
%% Get autospectra (for frequency finding)
AS{1} = squeeze(permMod{1}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS{2} = squeeze(permMod{2}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
% AS{3} = squeeze(permMod{3}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
clear('permMod','permModHD','permModSTR')

%% CROSS OVER
R.BBA_path = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\BetaBursts';
addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\statistics');
addpath('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\Plotting');
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\superbar'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\PEB_DFA');
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MLDFA');
load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_t5000_Sims.mat'])

% Compute Bursts
R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condname = {'Full','No HD','No M-STR'};
R.condcmap = linspecer(3);
R.cohband = 2;
BB = compute_BetaBursts_Simulated(R);

% Plot Burst Properties
subplot(1,2,1); set(gca, 'XScale', 'log')
plotBurstDurationHistogram(R,BB)
subplot(1,2,2); set(gca, 'XScale', 'log')
plotBurstAmplitudeHistogram(R,BB)

%% Do Beta Locked Analysis
addpath(R.BBA_path)

PLVcmp = 0;
close all
F(1) = figure; F(2) = figure;
BA_analysis(R,1,BB,xsimMod{1}{1},AS{1},PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model')

BA_analysis(R,2,BB,xsimMod{2}{1},AS{2},PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

% BA_analysis(R,3,BB,xsimMod{3}{1},AS{3},PLVcmp,F);
% figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

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

% BA_analysis(R,3,BB,xsimMod{3}{1},AS{3},PLVcmp,F);
% figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])

