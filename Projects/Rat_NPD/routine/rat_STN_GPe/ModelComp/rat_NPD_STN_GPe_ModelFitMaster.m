%STN GPE MOD FIT MASTER
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%  delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])

% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
closeMessageBoxes


rng(7564332)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;

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
        
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(modID)]);
        [R,p,m] = modelspec(R);
        pause(5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_STN_GPe_ModComp_M%.0f',modID); % 'All Cross'
        R.SimAn.rep = 512;
        R = setSimTime(R,48);
        R.Bcond = 0;
        [p] = SimAn_ABC_060219(R,p,m);
        closeMessageBoxes()
    end
end
