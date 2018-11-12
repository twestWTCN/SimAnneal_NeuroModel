% PLOT MODEL RELATIVE PROBABILITIES
clear; close all
%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Get Data Features
R = prepareRatData_InDirect_Group_NPD(R);

R.modcomp.modEvi.eps = -0.5;
R.modcomp.modN = 18;
R.modcompplot.NPDsel = [6 12 15];
plotModComp_091118(R)