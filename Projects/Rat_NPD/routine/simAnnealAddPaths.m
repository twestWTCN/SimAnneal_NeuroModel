function simAnnealAddPaths()

pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi('C:\spm12', pathCell));

if ~onPath; addpath('C:\spm12'); spm eeg; close all; end

if strmatch(getenv('computername'),'SFLAP-2')
    usname = 'Tim';
else
    usname = 'twest';
end

% addpath(['C:\Users\tim\Documents\Work\MATLAB ADDONS\parforprogmon')
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\aboxplot'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\ATvDFA-package'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\bplot\'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\Circular_Statistics_Toolbox'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\DrosteEffect-BrewerMap-221b913'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\export_fig'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\FMINSEARCHBND'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\linspecer'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\MEG_STN_Project'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\sigstar-master'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\sort_nat'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\SplitVec'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\TWtools'])
addpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\violin'])
addpath(['C:\Users\tim\Documents\Work\MATLAB ADDONS\violin'])
addpath(genpath(['C:\Users\' usname '\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD']))
addpath(genpath(['C:\spm12\toolbox\xjview96\xjview']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\boundedline-pkg']))
addpath(genpath(['C:\Users\' usname '\Documents\Work\MATLAB ADDONS\DylanMuir-ParforProgMon-9a1c257']))

