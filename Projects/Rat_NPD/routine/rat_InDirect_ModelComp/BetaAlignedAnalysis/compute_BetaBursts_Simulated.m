function [R,BB] = compute_BetaBursts_Simulated(R,xsimMod)
for cond = 1:size(R.condname,2)
    % Setup Data Structure
    if cond == 1
        vc_clean.trial{1} = xsimMod{1}{1}{1}{1}([1:4],:);
    elseif cond > 1
        vc_clean.trial{1} = xsimMod{cond}{1}([1:4],:);
    end
    vc_clean.fsample = 1/R.IntP.dt;
    vc_clean.time{1} = linspace(0,size(vc_clean.trial{1},2)/2000,size(vc_clean.trial{1},2));
    vc_clean.label = R.chsim_name;
    
    % Resample to workable resolution
    cfg = [];
    cfg.resamplefs = 200;
    vc_clean = ft_resampledata(cfg,vc_clean);
    if isequal(R.condname{cond},'Empirical')
        BB.TEmp = vc_clean.time{1};
    end
    % Set up some parameters
    R.bandinits = {'\alpha','\beta_1','\beta_2'};
    R.cohband = 2;
    R.BB.PLmeth = 'PPC';
    R.BB.decompmeth.type = 'wavelet';
    %     R.BB.decompmeth.filter.bwid = 2.5; % filter bandwidth
    R.BB.SW.winsize = 0.25;
    R.BB.SW.winover = 0.90;
    BB.powfrq = 21;
    BB.cohfrq = 21;
    BB.fsamp = vc_clean.fsample;
    R.BB.pairInd = [1 4]; % first is reference, second is main channel (e.g. M2 and STN)
    % Now do the decomposition
    BB = compute_SpectralDecomposition(R,BB,vc_clean,cond,0);
end

% Setup Time Vectors
BB.fsamp_sw = 1/(BB.SWTvec{1}(5)-BB.SWTvec{1}(4));
R.fsamp = BB.fsamp;
BB.T = linspace(0,length([BB.A{1}])/BB.fsamp,length([BB.A{1}]));
BB.TSw = linspace(0,length([BB.PLV{1}])/BB.fsamp_sw,length([BB.PLV{1}]));

disp('Assumes condition 6 is empirical!!')
BB.TSwEmp = linspace(0,length([BB.PLV{6}])/BB.fsamp_sw,length([BB.PLV{6}]));
% % Switch for Surrogates
surflag = 0; plotop = 1;

save([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'],'BB','R')
