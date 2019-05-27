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
R = prepareRatData_Demo_NPD(R);
WML = [];
try
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    disp('Making Mod List!!')
end

for modID = 1
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        
        %% Prepare Model
        modelspec = eval(['@MS_' R.out.tag '_model' num2str(modID)]);
        [R p m uc] = modelspec(R); % M! intrinsics shrunk"
%         pause(5)
        R.out.dag = sprintf([R.out.tag '_Model_ %.0f'],modID); % 'All Cross'
        R.Bcond = 0;
        %% Run ABC Optimization
        [p] = SimAn_ABC_220219b(R,p,m);
        closeMessageBoxes
    end
end
