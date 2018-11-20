% PLOT MODEL RELATIVE PROBABILITIES
clear; close all
%% Add Paths
% simAnnealAddPaths()
% ADAPT SO THAT .EPS is some prctile of all models - removes arbitrary
% selection and makes comparison a bit more understandable
%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;

%% Get Data Features
R = prepareRatData_STN_GPe_NPD(R);

