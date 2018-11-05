% Compute Beta Burst Aligned Analysis
% clear; close all
% simAnnealAddPaths()
R.rootn = 'C:\Users\Tim\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'InDrt_ModComp';
daglist = {'NPD_InDrt_ModComp_M1','NPD_InDrt_ModComp_M2','NPD_InDrt_ModComp_M3',...
    'NPD_InDrt_ModComp_M4', 'NPD_InDrt_ModComp_M5', 'NPD_InDrt_ModComp_M6',...
    'NPD_InDrt_ModComp_M7', 'NPD_InDrt_ModComp_M8'};

%     load([R.rootn 'outputs\' R.out.tag '\modeProbs.mat'])
% permMod = varo;
mnum = 6; % best model

% Load Model
load([R.rootn 'outputs R.out.tag '\' daglist{mnum} '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
m = varo;
% load modelfit
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelfit_' R.out.tag '_' daglist{mnum} '.mat'])
A = varo;
p = A.BPfit;
% Load Options
load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\R_' R.out.tag '_' daglist{mnum} '.mat'])
R = varo;
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
R = setSimTime(R,128);
[permMod xsimMod] = modelProbs(m.x,m,p,R);

close all
figure
plot(R.IntP.tvec_obs,xsimMod{1}{1}')
legend(R.chsim_name)

BA_analysis(R,xsimMod{1},permMod)


