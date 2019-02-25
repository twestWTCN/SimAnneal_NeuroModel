%STN GPE MOD FIT MASTER
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%  delete([R.rootn 'outputs\' R.out.tag '\MultiStartList.mat'])

% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
closeMessageBoxes
rng(7657456)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;

%% Prepare the data
R = prepareRatData_STN_GPe_NPD(R);

try
    load([R.rootn 'outputs\' R.out.tag '\MultiStartList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
    disp('Making Mod List!!')
end
N = 10; % number of starts
%% Prepare Model
for multiStart = N+1
    load([R.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
    if ~any(intersect(WML,multiStart))
        WML = [WML multiStart];
        save([R.rootn 'outputs\' R.out.tag '\MultiStartList'],'WML')
        disp('Writing to Mod List!!')
        fprintf('Now Fitting Mulitstart %.0f',multiStart)
        f = msgbox(sprintf('Fitting Multistart %.0f',multiStart));
        if multiStart == N+1; R = prepareRatData_STN_GPe_NPD(R,1,1); end
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
        [R,p,m] = modelspec(R);
        pause(5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
        % Delete Previous Saves
        delete([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\*'])
        R.SimAn.rep = 300;
        R = setSimTime(R,36);
        R.Bcond = 0;
        R.plot.flag = 1;
        SimAn_ABC_220219(R,p,m);
        closeMessageBoxes
    end
end
