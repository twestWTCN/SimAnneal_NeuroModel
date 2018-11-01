% COMPUTE MODEL RELATIVE PROBABILITIES
% clear; close all
simAnnealAddPaths()
R.rootn = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\Rat_NPD\';
R.out.tag = 'InDrt_ModComp';
daglist = {'NPD_InDrt_ModComp_M1','NPD_InDrt_ModComp_M2','NPD_InDrt_ModComp_M3',...
    'NPD_InDrt_ModComp_M4', 'NPD_InDrt_ModComp_M5', 'NPD_InDrt_ModComp_M6',...
    'NPD_InDrt_ModComp_M7', 'NPD_InDrt_ModComp_M8'};
% If concatanating to previously computed model comp structure
% %     load([R.rootn 'outputs\' R.out.tag '\modeProbs_v2.mat'])
% % permMod = varo;

for mnum =1:numel(daglist)
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
    R.Mfit = A;
    a = eval(['@MS_rat_InDirect_ModelComp_Model' num2str(mnum)]);
    [dum prior] = a(R);
    R.Mfit.prior = prior;
    % load parbank?
    load([R.rootn 'outputs\' R.out.tag '\' daglist{mnum} '\parBank_' R.out.tag '_' daglist{mnum} '.mat'])
    parBank =  varo;
    R = setSimTime(R,32);
    
 R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
R.analysis.BAA = 0; % Turn off BAA flag (time-locked analysis)
    parOptBank = parBank(1:end-1,1:R.SimAn.minRank);
%     if size(parOptBank,2)>2^10
%         parOptBank = parOptBank(:,1:2^10);
%     end
    
    if  size(parOptBank,2)>1
        R.parOptBank = parOptBank;
        
        R.obs.gainmeth = R.obs.gainmeth(1);
        figure(mnum);
        R.analysis.modEvi.N = 1000;
        permMod{mnum} = modelProbs(m.x,m,p,R);
    else
        permMod{mnum} = [];
    end
    saveMkPath([R.rootn 'outputs\' R.out.tag '\modeProbs_v2.mat'],permMod)
end