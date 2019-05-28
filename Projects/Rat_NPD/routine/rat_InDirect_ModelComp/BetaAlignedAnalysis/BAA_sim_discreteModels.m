function [R] = BAA_sim_discreteModels(R,modID,simtime)
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
feat{1} = permMod{1}.feat_rep;
% Now Modify
% Model 2 (weakened hyperdirect)
Pbase = permMod{1}.par_rep{1};
Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*0.01); % reduce connection by 90%
[r2mean,pnew,feat_sim,xsimMod{2}] = computeSimData(R,m,[],Pbase,0);
feat{2} = feat_sim;
% Model 3 (strengthened hyperdirect)
Pbase = permMod{1}.par_rep{1};
Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*3); % strengthen connection by 90%
[r2mean,pnew,feat_sim,xsimMod{3}] = computeSimData(R,m,[],Pbase,0);
feat{3} = feat_sim;

% Model 4 (weakened striato-pallidal)
Pbase = permMod{1}.par_rep{1};
Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(3,4))*0.01); % reduce connection by 95%
[r2mean,pnew,feat_sim,xsimMod{4}] = computeSimData(R,m,[],Pbase,0);
feat{4} = feat_sim;

% Model 5 (weakened striato-pallidal)
Pbase = permMod{1}.par_rep{1};
Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(3,4))*1.25); % strengthen connection by 90%
[r2mean,pnew,feat_sim,xsimMod{5}] = computeSimData(R,m,[],Pbase,0);
feat{5} = feat_sim;

%% Now Add the Empirical Data
close all
load([R.rootn 'data\Storage\L6_lesion_rat_020317.mat'])
cfg = [];
cfg.resamplefs = 1/R.IntP.dt;
FTdata = ft_resampledata(cfg,FTdata);

chind = [1 12 6 16]; % M2 STR GPe STN (order of simulated data)
xsimMod{6}{1} = zscore(FTdata.trial{1}(chind,:),0,2);
fsamp = FTdata.fsample;
R = prepareRatData_InDirect_Group_NPD(R)
feat{6} = R.data.feat_emp;
R.plot.outFeatFx(feat{1},feat(2:3),R.frqz,R,1,[])

R = simannealsetup_InDirect_ModelComp();
mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],'xsimMod','feat')
