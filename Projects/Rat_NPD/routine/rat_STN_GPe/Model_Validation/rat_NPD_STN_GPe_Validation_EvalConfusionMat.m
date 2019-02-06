%STN GPE MOD FIT MASTER
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%   delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])

% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));

rng(5353252)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
load([R.rootn 'outputs\' R.out.tag '\ConfData'],'RSimData')
confmatlist = combvec(1:3,1:3);
R.tmp.confmat = confmatlist;
R.comptype = 2;
R.data = RSimData.data;
%% Do the model probability computations
modelCompMaster(R,4:9,1:3)
% modelCompMaster(R,2,[1 3])

%% Plot the modComp results
R.modcomp.modN = 1:9;
R.modcompplot.NPDsel = [1:9];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
plotModComp_091118(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
