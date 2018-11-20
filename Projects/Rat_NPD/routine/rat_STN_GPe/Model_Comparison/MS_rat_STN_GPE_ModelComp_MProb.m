% COMPUTE MODEL RELATIVE PROBABILITIES
clear; close all
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));

%% Add Paths
% simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;

