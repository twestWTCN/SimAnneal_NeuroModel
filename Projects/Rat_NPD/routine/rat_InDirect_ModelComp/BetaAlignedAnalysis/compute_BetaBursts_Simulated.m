function [R,BB] = compute_BetaBursts_Simulated(R,xsimMod)
for cond = 1:size(R.condname,2)
    % Setup Data Structure
    if cond == 1
        vc_clean.trial{1} = xsimMod{1}{1}{1}([1:4],:);
    elseif cond == 2
        vc_clean.trial{1} = xsimMod{2}{1}{1}([1:4],:);
    elseif cond == 3
        vc_clean.trial{1} = xsimMod{3}{1}{1}([1:4],:);
    end
    vc_clean.fsample = 2000;
    vc_clean.time{1} = linspace(0,size(xsimMod{1}{1}{1},2)/2000,size(xsimMod{1}{1}{1},2));
    vc_clean.label = R.chsim_name;
    
    % Resample to workable resolution
    cfg = [];
    cfg.resamplefs = 200;
    vc_clean = ft_resampledata(cfg,vc_clean);
    
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

% % Switch for Surrogates
surflag = 0; plotop = 1;

% Find the amplitude threshold from the prctile of concatanated data
if surflag == 0
    BB.epsAmp = prctile(BB.A{1}(R.BB.pairInd(2),:),80,2);
    BB.epsPLV = prctile(BB.PLV{1}(1,:),75,2);
else
    BORG = load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_org'],'BB');
    BB.epsAmp = BORG.BB.epsAmp;
    BB.epsPLV = BORG.BB.epsPLV ;
    clear BORG
end
% Plot Bursting
if plotop == 1
    figure(1)
    plotExampleBurstPLV(R,BB)
    set(gcf,'Position',[680  358  1048  622])
end

save([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'],'BB','R')
