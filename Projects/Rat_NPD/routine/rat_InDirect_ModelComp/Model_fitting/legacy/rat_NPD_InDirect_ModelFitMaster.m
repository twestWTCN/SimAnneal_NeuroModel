% MASTER FUNCTION FOR MODEL FITTING
% Models are specified in ModelSpecs of project folder. Will batch fit
% models to the data specified in fx prepareRatData. WML stores batch
% progress and can be accessed by multiple MATLAB realizations.
%%%%%%%%%%%%%%%%%%%%%%%%
% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));

% Add relevant paths for toolboxes
% simAnnealAddPaths()
rng(123213)

%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
% IF FRESH START
% delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])


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

for modID = [1:18]
    load([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingModList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Model %.0f',modID)
        f = msgbox(sprintf('Fitting Model %.0f',modID));
        
        %% Prepare Model
        modelspec = eval(['@MS_rat_' R.out.tag '_Model' num2str(modID)]);
        [R p m uc] = modelspec(R);
        pause(5)
        R.out.dag = sprintf('NPD_InDrt_ModComp_M%.0f',modID); % 'All Cross'
        
        %% Run ABC Optimization
        R = setSimTime(R,32);
        R.Bcond = 0;
        parBank = [];
        R.SimAn.rep = 512; %448
        [p] = SimAn_ABC_211218(m.x,uc,p,m,R,parBank);
    end
end
