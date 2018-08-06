clear; close all
R = simannealsetup_STN_GPe_M2_ModelComp();
[dum prior]  = MS_rat_STN_GPe_M2_ModelComp_Model6_wMd2Priors(R);

load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\ModComp\Mod6\modelfit_ModComp_NPD_ModComp_M6.mat')
R =varo;
p = R.Mfit.BPfit;
R.Mfit.prior = prior;
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\ModComp\Mod6\modelspec_ModComp_NPD_ModComp_M6.mat')
m = varo;
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\outputs\ModComp\Mod6\parBank_ModComp_NPD_ModComp_M6.mat')
parBank = varo;
R.analysis.modEvi.eps = -0.35;
R.analysis.modEvi.N = 100;


[KL DKL] = KLDiv(R,p,m,parBank);
