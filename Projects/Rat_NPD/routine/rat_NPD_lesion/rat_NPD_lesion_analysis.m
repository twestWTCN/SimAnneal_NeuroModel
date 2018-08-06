clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\boundedline-pkg'))
R = simannealsetup_NPD_lesion_rat;
R.out.dag = '0520818_lesion_cross'; % cross only

load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
% load modelfit
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
p = R.Mfit.BPfit;
% load parbank?
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;


R = setSimTime(R,128);
R.analysis.modEvi.N = 100;
R.analysis.modEvi.eps = prctile(parBank(end,:),90);

ABC_Gen_Analysis(R,p,m,parBank)
