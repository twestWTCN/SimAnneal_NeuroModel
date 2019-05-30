function burstAnalyAddPaths()
switch getenv('computername')
    case 'SFLAP-2'
        usname = 'Tim'; gitpath = '\Documents\Work\GIT'; madpath = 'MATLAB_ADDONS';
        spmpath = 'C:\Users\Tim\Documents\spm12';
    case 'FREE'
        usname = 'twest'; gitpath = '\Documents\Work\GitHub'; madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\spm12';
    case 'DESKTOP-94CEG1L'
        usname = 'timot'; gitpath =  '\Documents\GitHub'; madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
end

pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end
addpath(['C:\Users\' usname '\Documents\' madpath '\ParforProgMon'])
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
addpath(genpath(['C:\Users\' usname '\Documents\' madpath '\ParforProgMon']))
addpath(['C:\Users\' usname '\Documents\' madpath '\sigstar-master'])
addpath(['C:\Users\' usname '\Documents\' madpath '\sort_nat'])
addpath(['C:\Users\' usname '\Documents\' madpath '\SplitVec'])
addpath(['C:\Users\' usname '\Documents\' madpath '\TWtools'])
addpath(['C:\Users\' usname '\Documents\' madpath '\violin'])
addpath(genpath(['C:\spm12\toolbox\xjview96\xjview']))
addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\sim_machinery']))
addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_InDirect_ModelComp']))
addpath(genpath(['C:\Users\' usname '\Documents\' madpath '\boundedline-pkg']))
addpath(genpath(['C:\Users\' usname '\' gitpath '\BrewerMap']))
addpath(genpath(['C:\Users\' usname '\' gitpath '\BurstToolbox']))
addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\data']));
addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\model_fx']));
addpath(genpath(['C:\Users\' usname '\' gitpath '\SimAnneal_NeuroModel\Projects\Rat_NPD\priors']));