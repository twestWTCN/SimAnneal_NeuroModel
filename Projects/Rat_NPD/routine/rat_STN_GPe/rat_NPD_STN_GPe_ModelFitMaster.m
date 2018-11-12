%STN GPE MOD FIT MASTER
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
% delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])
%

% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));

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
rng(4342131)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
R.objfx.specspec = 'cross'; 
%% Prepare the data
R = prepareRatData_STN_GPe_NPD(R);

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
for modID = 1:3
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        
        modelspec = eval(['@MS_rat_STN_GPe_ModelComp_Model' num2str(modID)]);
        [R p m uc] = modelspec(R);
        pause(5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_STN_GPe_ModComp_M%.0f',modID); % 'All Cross'
        
        R = setSimTime(R,32);
        R.Bcond = 0;
        parBank = [];
        [p] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);
    end
end
