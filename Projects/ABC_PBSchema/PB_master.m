clear; close all
% MASTER FOR PB SCHEMA
R = PB_scheme_setup();
R = addPaths_PBSchema(R);

% load('C:\Users\Tim\Dropbox\Flavie_project\ABC_PBSchema\Data\schema2data')
% R.data.feat_emp = feat_emp;

[R] = getFlaviesData(R);

PB_PlotOutputs(R,{R.data.feat_emp},[],NaN);
[R p m uc] = setup_model_testing(R);

% Precompute Beta Table
% precomputeBetaSig(R,m)
% precomputeTMS(R,m)
% 
R.plot.flag = 1;
R.plot.save = 0;

[p] = SimAn_ABC_220219(R,p,m);

% LOOK AT REDUCING THE COV OF THE REAL DATA - THE EFFECT IS THERE BUT
% SHROUDED IN NOSIE!