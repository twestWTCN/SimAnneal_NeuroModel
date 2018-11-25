% COMPUTE MODEL RELATIVE PROBABILITIES
% clear; close all
%% Add Paths
simAnnealAddPaths()
%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;
modelCompMaster(R,WML)

%% Setup for parallelisation (multiple MATLAB sessions)
try
    load([R.rootn 'outputs\' R.out.tag '\WorkingPermModList'])
    disp('Loaded Perm Mod List!!')
    % If concatanating to previously computed model comp structure
catch
    WML = [];
    mkdir([R.rootn 'outputs\' R.out.tag ]);
    save([R.rootn 'outputs\' R.out.tag '\WorkingPermModList'],'WML')
    disp('Making Perm Mod List!!')
end

%% Main Loop
for modID =1:18
    load([R.rootn 'outputs\' R.out.tag '\WorkingPermModList'],'WML')
    permMod = [];
    if ~any(intersect(WML,modID))
        WML = [WML modID];
        save([R.rootn 'outputs\' R.out.tag '\WorkingPermModList'],'WML')
        disp('Writing to PermMod List!!')
        fprintf('Now Computing Probabilities for Model %.0f',modID)
        f = msgbox(sprintf('Probabilities for Model %.0f',modID));
       
        
        dagname = sprintf('NPD_InDrt_ModComp_M%.0f',modID);
        % Load Model
        load([R.rootn 'outputs\' R.out.tag '\' dagname '\modelspec_' R.out.tag '_' dagname '.mat'])
        m = varo;
        % load modelfit
        load([R.rootn 'outputs\' R.out.tag '\' dagname '\modelfit_' R.out.tag '_' dagname '.mat'])
        A = varo;
        p = A.BPfit;
        % Load Options
        load([R.rootn 'outputs\' R.out.tag '\' dagname '\R_' R.out.tag '_' dagname '.mat'])
        R = varo;
        R.Mfit = A;
        a = eval(['@MS_rat_InDirect_ModelComp_Model' num2str(modID)]);
        [dum prior] = a(R);
        R.Mfit.prior = prior;
        % load parbank?
        load([R.rootn 'outputs\' R.out.tag '\' dagname '\parBank_' R.out.tag '_' dagname '.mat'])
        parBank =  varo;
        R = setSimTime(R,32);
        
        R.analysis.modEvi.eps = parBank(end,R.SimAn.minRank);
        R.analysis.BAA = 0; % Turn off BAA flag (time-locked analysis)
        parOptBank = parBank(1:end-1,1:R.SimAn.minRank);
        %     if size(parOptBank,2)>2^10
        %         parOptBank = parOptBank(:,1:2^10);
        %     end
        
        if  size(parOptBank,2)>1
            R.parOptBank = parOptBank;
            R.obs.gainmeth = R.obs.gainmeth(1);
            figure(modID);
            R.analysis.modEvi.N = 1000;
            permMod = modelProbs(m.x,m,p,R);
        else
            permMod = [];
        end
        saveMkPath([R.rootn 'outputs\' R.out.tag '\' dagname '\modeProbs_' R.out.tag '_' dagname '.mat'],permMod)
        pause(10)
        close all
    end
end