clear ; close all
R = simannealsetup_NPD_300817_OFF;
addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\boundedline-pkg'))
d = 'trialanalysisOFF';
% Need m,p, parBank
% load('emerg_save.mat','p')

% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\M_OFF_FRESH.mat')
% 
R.out.tag = 'NPD_Final_JNPPaper';
R.out.dag = '180718_COND1_auto';
R = simannealsetup_NPD_Can_060718()
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913')
% load('C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\FULL_OFF_CROSS.mat')
% p = bmod;

 R.out.dag = '180718_COND1_full_run2'; % cross only
R.plot.save = 'True';
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
xobs1 = varo;
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;
p = spm_unvec(parBank(1:end-1,1),xobs1.Mfit.Pfit);
load(['C:\Users\twest\Documents\Work\GitHub\Phase_Sync_and_Stability\Models\NMM\ABC_Fitted\Full_network\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])

R = setSimTime(R,128);
ABC_Gen_Analysis(R,p,m,parBank)
