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
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        [R p m uc] = MS_rat_InDirect_ModelComp_Model1(R);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_InDrt_ModComp_M%.0f',modID); % 'All Cross'
        
        R = setSimTime(R,32);
        R.Bcond = 0;
        parBank = [];
        [p] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);
    end
end
