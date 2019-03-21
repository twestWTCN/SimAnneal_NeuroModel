clear; close all
% MASTER FOR PB SCHEMA
R = PB_scheme_setup();
% R = addPaths_PBSchema(R)

% load('C:\Users\Tim\Dropbox\Flavie_project\ABC_PBSchema\Data\schema2data')
% R.data.feat_emp = feat_emp;

[R] = getFlaviesData(R);

PB_PlotOutputs(R,{R.data.feat_emp},[],NaN);
[R p m uc] = setup_model_testing(R);

% Precompute Beta Table
precomputeBetaSig(R,m)
precomputeTMS(R,m)

R.plot.flag = 1;
R.plot.save = 0;

[p] = SimAn_ABC_220219(R,p,m);

varo = [];
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'])
pfit = spm_unvec(varo(1:end-1,1),p);


u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
[r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,pfit,0,0);

