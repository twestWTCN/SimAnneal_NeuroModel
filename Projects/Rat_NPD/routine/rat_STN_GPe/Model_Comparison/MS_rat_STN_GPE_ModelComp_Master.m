% COMPUTE MODEL RELATIVE PROBABILITIES AND PLOT RESULTS
clear; close all
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));


%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
% Get empirical data
R = prepareRatData_STN_GPe_NPD(R,0); 

% Do the model probability computations
% modelCompMaster(R)

% Plot the modComp results
R.modcomp.modN = 3;
R.modcompplot.NPDsel = [1:3];
R.plot.confint = 'yes';
cmap = linspecer(R.modcomp.modN);
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
