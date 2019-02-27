% COMPUTE MODEL RELATIVE PROBABILITIES AND PLOT RESULTS
clear; close all
closeMessageBoxes

%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
% Get empirical data
R = prepareRatData_InDirect_Group_NPD(R); 

%% Do the model probability computations
R.comptype = 1;
modelCompMaster(R,8:9) %,[1:8 10:12]

%% Plot the modComp results
R.modcomp.modN = [1:12];
R.modcompplot.NPDsel = [4 6 10]; %[6 9 10];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
cmap = cmap(end:-1:1,:);
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
