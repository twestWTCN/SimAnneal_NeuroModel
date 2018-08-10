% simAnnealAddPaths()
close all

daglist = {'NPD_ModComp_M1','NPD_ModComp_M2','NPD_ModComp_M3',...
    'NPD_ModComp_M4','NPD_ModComp_M5','NPD_ModComp_M6','NPD_ModComp_M7'};
R.rootn = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'ModComp';

load([R.rootn 'outputs\' R.out.tag '\' daglist{1} '\modelfit_' R.out.tag '_' daglist{1} '.mat'])
R = varo;
load([R.rootn 'outputs\' R.out.tag '\modeProbs.mat'])
permMod = varo;

 R.analysis.modEvi.eps = -0.2;
plotModComp(R,permMod)