function simAnnealAddPaths()

pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi('C:\spm12', pathCell));


if strmatch(getenv('computername'),'SFLAP-2')
    usname = 'Tim'; gitname = 'GIT'; madpath = 'MATLAB_ADDONS';
        if ~onPath; addpath('C:\Users\Tim\Documents\spm12'); spm eeg; close all; end
else
    usname = 'twest'; gitname = 'GitHub'; madpath = 'Work\MATLAB ADDONS';
    if ~onPath; addpath('C:\spm12'); spm eeg; close all; end
end

% addpath(['C:\Users\tim\Documents\Work\MATLAB ADDONS\parforprogmon')
addpath(['C:\Users\' usname '\Documents\' madpath '\aboxplot'])
addpath(['C:\Users\' usname '\Documents\' madpath '\ATvDFA-package'])
addpath(['C:\Users\' usname '\Documents\' madpath '\bplot\'])
addpath(['C:\Users\' usname '\Documents\' madpath '\Circular_Statistics_Toolbox'])
addpath(['C:\Users\' usname '\Documents\' madpath '\DrosteEffect-BrewerMap-221b913'])
addpath(['C:\Users\' usname '\Documents\' madpath '\export_fig'])
addpath(['C:\Users\' usname '\Documents\' madpath '\FMINSEARCHBND'])
addpath(['C:\Users\' usname '\Documents\' madpath '\linspecer'])
addpath(['C:\Users\' usname '\Documents\' madpath '\MEG_STN_Project'])
addpath(['C:\Users\' usname '\Documents\' madpath '\Neurospec\neurospec21'])
addpath(['C:\Users\' usname '\Documents\' madpath '\sigstar-master'])
addpath(['C:\Users\' usname '\Documents\' madpath '\sort_nat'])
addpath(['C:\Users\' usname '\Documents\' madpath '\SplitVec'])
addpath(['C:\Users\' usname '\Documents\' madpath '\TWtools'])
addpath(['C:\Users\' usname '\Documents\' madpath '\violin'])
addpath(genpath(['C:\spm12\toolbox\xjview96\xjview']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\' gitname '\SimAnneal_NeuroModel\sim_machinery']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\' gitname '\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_InDirect_ModelComp']))
addpath(genpath(['C:\Users\' usname '\Documents\' madpath '\boundedline-pkg']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\' gitname '\BrewerMap']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\' gitname '\BurstToolbox']))
