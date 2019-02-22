% COMPUTE MODEL RELATIVE PROBABILITIES AND PLOT RESULTS
clear; close all
        closeMessageBoxes

%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
% Get empirical data
R = prepareRatData_STN_GPe_NPD(R,0); 

%% Do the model probability computations
R.comptype = 1;
modelCompMaster(R,1:3,[])

%% Plot the modComp results
R.modcomp.modN = 1:3;
R.modcompplot.NPDsel = [1:3];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])

