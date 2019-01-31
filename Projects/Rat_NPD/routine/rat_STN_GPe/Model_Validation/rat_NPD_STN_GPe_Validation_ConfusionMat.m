%STN GPE MOD FIT MASTER
%%%%%%%%%%%%%%%%%%%%%%%%
% IF FRESH!
%   delete([R.rootn 'outputs\' R.out.tag '\WorkingModList.mat'])

% simAnnealAddPaths()
clear ; close all

% Close all msgboxes
handles = allchild(0);
tags = get(handles,'Tag');
isMsg = strncmp(tags,'Msgbox_',7); % all message boxes have the tags in the format of Msgbox_*
delete(handles(isMsg));

rng(5353252)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;

for modID = 1:3
    Rt = []; % Temp R Struc
    %% First Simulate the Data from the Empircally Fitted Models
    % Recover Fitted Parameters
    dagname = sprintf([R.out.tag '_M%.0f'],modID);
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modelfit_' R.out.tag '_NPD_' dagname '.mat'])
    Mfit = varo; %i.e. permMod
    p = Mfit.Pfit;
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modelspec_' R.out.tag '_NPD_' dagname '.mat'])
    m = varo; %i.e. permMod
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\R_' R.out.tag '_NPD_' dagname '.mat'])
    Rt = varo; %i.e. permMod
    Rt.root = R.rootn; %Overwrite root with current
    
    % Retrieve structure map
    pInd = parOptInds_110817(Rt,p,m.m); % in structure form
    pIndMap = spm_vec(pInd); % in flat form
    
    % Draw copula and average
    par = postDrawCopula(Mfit,p,pIndMap,1000);
    avPar = averageCell(par);
    
    % Now Simulate Fitted Model
    pnew = avPar;
    [r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData(Rt,m,[],pnew,128,1);
    
    Rt.data.feat_emp = feat_sim;
    % squeeze(meannpd_data(1,1,1,1,:))
    Rt.data.feat_xscale = R.frqz;
    RSimData(modID) = Rt;
end
save([R.rootn 'outputs\' R.out.tag '\ConfData'],'RSimData')

load([R.rootn 'outputs\' R.out.tag '\ConfData'],'RSimData')
% Make matrix of combinations for confusion matrix
confmatlist = combvec(1:3,1:3);

% Create List (for parallelization across multiple MATLAB instances)
delete([R.rootn 'outputs\' R.out.tag '\ConfWorkList'])

try
    load([R.rootn 'outputs\' R.out.tag '\ConfWorkList'])
    disp('Loaded Mod List!!')
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\ConfWorkList'],'WML')
    disp('Making Mod List!!')
end

%% Prepare Model
for i = 1:size(confmatlist,2)
    load([R.rootn 'outputs\' R.out.tag '\ConfWorkList'],'WML')
    if ~any(intersect(WML,i))
        WML = [WML i];
        save([R.rootn 'outputs\' R.out.tag '\ConfWorkList'],'WML')
        disp('Writing to Mod List!!')
        
        SimData = confmatlist(1,i)
        SimMod = confmatlist(1,i)
        
        
        fprintf('Fitting Model %.0f to data %.0f',SimMod,SimData)
        f = msgbox(sprintf('Fitting Model %.0f to data %.0f',SimMod,SimData));
        
        modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(SimMod)]);
        [R p m uc] = modelspec(R);
        pause(5)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R.data = RSimData(SimData).data;
        R.out.dag = sprintf('NPD_STN_GPe_ConfMat_DataM%.0f_ParM%.0f',SimData,SimMod); % 'All Cross'
        R.SimAn.rep = 64;
        R = setSimTime(R,32);
        R.Bcond = 0;
        parBank = [];
        [p] = SimAn_ABC_211218(m.x,uc,p,m,R,parBank);
    end
end
