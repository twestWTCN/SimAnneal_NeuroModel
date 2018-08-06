% COMPUTE MODEL RELATIVE PROBABILITIES
close all
R.out.tag = 'ModComp';
daglist = {'NPD_ModComp_M1','NPD_ModComp_M2','NPD_ModComp_M3',...
    'NPD_ModComp_M4','NPD_ModComp_M5','NPD_ModComp_M6'};

for mnum = 1:6
    % Load Model
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
    m = varo;
    % load modelfit
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
    R = varo;
    p = R.Mfit.BPfit;
    % load parbank?
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' daglist{mnum} '.mat'])
    parBank =  varo;
    
    R.analysis.modEvi.eps = -0.35;
    parOptBank = parBank(1:end-1,parBank(end,:)>R.analysis.modEvi.eps);
    R.parOptBank = parOptBank;
    
    R.obs.gainmeth = R.obs.gainmeth(1);
    figure(2);
    R.analysis.modEvi.N = 100;
    permMod{mnum} = modelProbs(m.x,m,p,R);
    saveMkPath([R.rootn 'outputs\' R.out.tag '\modeProbs.mat'],permMod)
end