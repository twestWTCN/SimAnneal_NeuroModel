function compute_BetaBursts_Simulated(R)

load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\routine\rat_InDirect_ModelComp\BetaBurstAnalysis\Data\Indirect_Sims.mat')
R = makeHeader_SubCort_Cort_Networks();
vc_clean.trial{1} = xsimMod{1}{1}([1 4],:);
vc_clean.fsample = 2000;
vc_clean.time{1} = linspace(0,size(xsimMod{1}{1},2)/2000,size(xsimMod{1}{1},2));
vc_clean.label = {'MMC','STN'};

cfg = [];
cfg.resamplefs = 200;
vc_clean = ft_resampledata(cfg,vc_clean);
%%%
close all
plotop = 1;
R.BB.PLmeth = 'PPC';
R.bandname = {'Alpha','B1','B2'};
% surflag = 0;
% For surrplotting
% R.condcmap(3,:) = [0.7 0.7 0.7];
% R.condname{3} = 'surrogate';
R.bandef = [8 12;16 24;24 36];
BB.powfrq = 18;
BB.cohfrq = 18;
BB = compute_BetaBursts(R,BB,vc_clean,1,0);
% BB.APrc{cond}
BB.fsamp = vc_clean.fsample;

F = []; % Empty Figure Handle

BB.fsamp_sw = 1/(BB.SWTvec{1}(5)-BB.SWTvec{1}(4));
R.fsamp = BB.fsamp;
BB.T = linspace(0,length([BB.A{1:2}])/BB.fsamp,length([BB.A{1:2}]));
BB.TSw = linspace(0,length([BB.PLV{1:2}])/BB.fsamp_sw,length([BB.PLV{1:2}]));

% hack for surrogate plots inside single subj
%                                         surflag = 0;
% ThresholdAmplitude
if surflag == 0
    BB.epsAmp = prctile([BB.A{:}],65,2);
    BB.epsPLV = prctile([BB.PLV{:}],65,2);
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
%% Now Do Stats and Plot Single Example (if desired)
% Amplitude Only
%                     BB.range.Amp = -100:5:100;
%                     BB.range.segDur = 3:150:1800;
BB.range.Amp = linspace(0,50,12); %0:5:120; % Group: 0:3:120; single: 0:5:80
BB.range.segDur = linspace(0,1800,12); %0:100:1800; % Group: 0:100:1800; single: 25:150:1800
BB.range.AmpPrc = 0:5:100;

BB.plot.lims.PLV = [-100 100]; %[0 0.35];
BB.plot.lims.wPLV = [-15 10];
BB.plot.lims.Amp =  [0 40];
BB.plot.lims.wAmpPrc =  [-2 8];
BB.plot.lims.Dur = [0 2500];

if plotop == 1; F = figure(2); end
BB = compute_BetaBurstsStats(R,BB,F,plotop);
if plotop == 1; set(F(1),'Position',[680.0000   83.5000  907.5000  894.5000]); end

% Amplitude and Sync
BB.range.PLV = linspace(-100,100,8); %linspace(0,0.36,8); %0:0.1:1;
if plotop == 1; F(1) = figure(3); end
BB = computeBetaBurstPLVStats(R,BB,F,plotop);

% Relative Phase
if plotop == 1; F(1) = figure(3);set(F(1),'Position',[638 291 1201 684]); end
BB.range.RP = linspace(-pi,pi,8); %-pi:pi/3:pi
BB = computeBetaBurstRPStats(R,BB,F,plotop);

% Save to data structure
%     mkdir([R.datapathr R.subname{sub} '\ftdata\BetaBursts'])
%     save([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_' surrtag],'BB')

% ! shutdown /h