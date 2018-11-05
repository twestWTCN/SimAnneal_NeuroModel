% MODEL 1:
% SERIAL FLOW
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

try
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    disp('Making Mod List!!')
end
%% Prepare Model
for modID = 1:14
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        
        [R p m uc] = MS_rat_InDirect_ModelComp_Model1(R);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_InDrt_ModComp_M%.0f',modID); % 'All Cross'
        R.out.tag = 'InDrt_ModComp';
        R.plot.save = 'False';
        R.obs.logdetrend =1;
        R.obs.trans.norm = 1;
        
        R.SimAn.starttemp = 2;
        R.SimAn.rep =448; %256; %96; %512; % Repeats per temperature
        R.SimAn.Tm = 1;
        R.SimAn.jitter = 0.5;
        R.SimAn.searchN = 200;
        R = setSimTime(R,32);
        R.objfx.specspec = 'cross';
        R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A','.S','.D'}; %,'.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
        R.Bcond = 0;
        parBank = [];
        R.SimAn.copout = [2 3];
        [p] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);
    end
end
