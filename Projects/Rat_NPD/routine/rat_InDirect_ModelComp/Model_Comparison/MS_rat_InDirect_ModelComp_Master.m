% COMPUTE MODEL RELATIVE PROBABILITIES AND PLOT RESULTS
clear; close all
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));


%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
% Get empirical data
R = prepareRatData_InDirect_Group_NPD(R); 

%% Do the model probability computations
% modelCompMaster(R,[100 200],[1:18])

%% Plot the modComp results
R.modcomp.modN = [1:18 100 200];
R.modcompplot.NPDsel = [6 12 15];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
