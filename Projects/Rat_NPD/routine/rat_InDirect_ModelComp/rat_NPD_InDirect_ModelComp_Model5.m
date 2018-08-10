% MODEL 5:
% SERIAL FLOW with Hyperdirect + M2 Auto Priors
%%%%%%%%%%%%%%%%%%%%%%%%
% simAnnealAddPaths()
clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery'))
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools\')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot\')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
addpath('C:\spm12')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
rng(24312321)

%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

%% Prepare the data
R = prepareRatData_InDirect_Group_NPD(R);

%% Prepare Model
[R p m uc] = MS_rat_InDirect_ModelComp_Model5(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R.out.dag = 'NPD_InDrt_ModComp_M5'; % 'All Cross'
R.out.tag = 'InDrt_ModComp';
R.plot.save = 'True';
R.obs.trans.norm = 1;
R.obs.logdetrend =1;

load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
p = R.Mfit.BPfit;
R.SimAn.starttemp = 1.5;

R.SimAn.rep =256; %256; %96; %512; % Repeats per temperature
R.SimAn.Tm = 1;
R.SimAn.jitter = 1.25;
R.SimAn.searchN = 200;
R = setSimTime(R,18);
R.objfx.specspec = 'cross';
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A','.S','.D'}; %,'.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.Bcond = 0;
parBank = [];
R.SimAn.copout = [2 3];
[p] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);

