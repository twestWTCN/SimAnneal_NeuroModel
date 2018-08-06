% COMPUTE MODEL RELATIVE PROBABILITIES
% clear; close all
R.rootn = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'ModComp';
daglist = {'NPD_ModComp_M1','NPD_ModComp_M2','NPD_ModComp_M3',...
    'NPD_ModComp_M4','NPD_ModComp_M5','NPD_ModComp_M6','NPD_ModComp_M7'};
for mnum =7
    % Load Model
    load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
    m = varo;
    % load modelfit
    load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\modelfit_' R.out.tag '_' daglist{mnum} '.mat'])
    R = varo;
    p = R.Mfit.BPfit;
    % load parbank?
    load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\parBank_' R.out.tag '_' daglist{mnum} '.mat'])
    parBank =  varo;
    R = setSimTime(R,20);
    
    R.analysis.modEvi.eps = -0.2;
    parOptBank = parBank(1:end-1,parBank(end,:)>R.analysis.modEvi.eps);
    if size(parOptBank,2)>2^10
        parOptBank = parOptBank(:,1:2^10);
    end
    
    if  size(parOptBank,2)>1
        R.parOptBank = parOptBank;
        
        R.obs.gainmeth = R.obs.gainmeth(1);
        figure(mnum);
        R.analysis.modEvi.N = 200;
        permMod{mnum} = modelProbs(m.x,m,p,R);
    else
        permMod{mnum} = [];
    end
    saveMkPath([R.rootn 'outputs\' R.out.tag '\modeProbs.mat'],permMod)
end