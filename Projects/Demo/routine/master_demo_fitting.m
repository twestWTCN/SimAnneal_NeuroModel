% MASTER- DEMO FOR MODEL FITTING
% This function incorporates all of the nodes (including GPi/Thalamus)
clear ; close all

% Add relevant paths for toolboxes
demoAddPaths()
% Close all msgboxes
closeMessageBoxes

%% Set Routine Pars
R = setup_demo;
R = setSimTime(R,32);

% IF FRESH START
% delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])
%% Prepare the data
% R = prepareRatData_Demo_NPD(R);

WML = [];
save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')

try
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    disp('Making Mod List!!')
end

% Model List:
% Model 0: Fit just the intrinsic dynamics of cortex
% Model 1: Full Model using parameters fit from model 0;

for modID = 2
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    if ~any(intersect(WML,modID))
        if modID == 0
            R.data.chsel = 1;
            R = prepareRatData_Demo_NPD(R);
            R.chloc_name = {'MMC'};
            R.chsim_name = {'MMC'};
        elseif modID == 1
            R.data.chsel = [1 2 4]';
            R.obs.obsstates = [1 2 3 4];
            R = prepareRatData_Demo_NPD(R);
            R.chloc_name = {'MMC'  'GPE'    'STN'};
            R.chsim_name = {'MMC'  'STR'    'GPE'    'STN'};
        elseif modID == 2
            R.data.chsel = [2 4]';
            R.obs.obsstates = [1 2 3];
            R = prepareRatData_DemoSynthPeak_NPD(R);
            R.chloc_name = {'GPE'    'STN'};
            R.chsim_name = {'GPE'    'STN'  'GPI'};
            
        end
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        
        %% Prepare Model
        modelspec = eval(['@MS_' R.out.tag '_M' num2str(modID)]);
        [R p m uc] = modelspec(R); % M! intrinsics shrunk"
        %         pause(5)
        R.out.dag = sprintf([R.out.tag '_M%.0f'],modID); % 'All Cross'
        R.Bcond = 0;
        %% Run ABC Optimization
        [p] = SimAn_ABC_220219b(R,p,m);
        closeMessageBoxes
    end
end
