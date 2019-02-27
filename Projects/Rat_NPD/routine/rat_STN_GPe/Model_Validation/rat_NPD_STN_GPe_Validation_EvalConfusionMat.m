%STN GPE CONFUSION MATRIX EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%
rat_NPD_STN_GPe_ModelFitMaster
% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
closeMessageBoxes
rng(5353252)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
load([R.rootn 'outputs\' R.out.tag '\ConfData'],'RSimData')
confmatlist = combvec(1:3,1:3);
R.tmp.confmat = confmatlist;
R.comptype = 2;
R.data = RSimData.data;
%% Do the model probability computations
modelCompMaster(R,1:9,[])
% modelCompMaster(R,2,[1 3])

%% Plot the modComp results
R.modcomp.modN = 1:9;
R.modcompplot.NPDsel = [1:9];
R.plot.confint = 'yes';
cmap = linspecer(numel(R.modcomp.modN));
plotModValidationEval(R,cmap)
figure(2)
subplot(3,1,1); ylim([-2 1])
