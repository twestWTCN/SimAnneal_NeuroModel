% Compute Beta Burst Aligned Analysis
% clear; close all
% simAnnealAddPaths()
% R.rootn = 'C:\Users\Tim\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.rootn = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\';

R.out.tag = 'InDrt_ModComp';
daglist = {'NPD_InDrt_ModComp_M1','NPD_InDrt_ModComp_M2','NPD_InDrt_ModComp_M3',...
    'NPD_InDrt_ModComp_M4', 'NPD_InDrt_ModComp_M5', 'NPD_InDrt_ModComp_M6',...
    'NPD_InDrt_ModComp_M7', 'NPD_InDrt_ModComp_M8'};

%     load([R.rootn 'outputs\' R.out.tag '\modeProbs.mat'])
% permMod = varo;
mnum = 6; % best model

% Load Model
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
m = varo;
% load modelfit
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelfit_' R.out.tag '_' daglist{mnum} '.mat'])
A = varo;
p = A.BPfit;
% Load Options
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\R_' R.out.tag '_' daglist{mnum} '.mat'])
R = varo;
R.rootn = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\';

R.Mfit = A;

a = eval(['@MS_rat_InDirect_ModelComp_Model' num2str(mnum)]);
[dum prior] = a(R);
R.Mfit.prior = prior;
% load parbank?
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\parBank_' R.out.tag '_' daglist{mnum} '.mat'])
parBank =  varo;

R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
parOptBank = parBank(1:end-1,1:R.SimAn.minRank);

R.parOptBank = parOptBank;
R.obs.gainmeth = R.obs.gainmeth(1);
figure;
R.analysis.modEvi.N = 1000;
R.analysis.BAA = 1;
R = setSimTime(R,1000);

% With Hyperdirect
p1 = p;
[permMod xsimMod] = modelProbs(m.x,m,p1,R);

% No Hyperdirect
p2 = p;
p2.A{1}(4,1) = -32;
[permModHD xsimModHD] = modelProbs(m.x,m,p2,R);

% No Hyperdirect
p3 = p;
p3.A{1}(2,1) = -32;
[permModSTR xsimModSTR] = modelProbs(m.x,m,p3,R);

close all
subplot(1,3,1)
BA_analysis(R,xsimMod{1},permMod)

subplot(1,3,2)
BA_analysis(R,xsimModHD{1},permModHD)

subplot(1,3,3)
BA_analysis(R,xsimModSTR{1},permModSTR)
