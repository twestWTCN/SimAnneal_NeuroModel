function BB = AmpDurStatistics(R,BB)
% Setup Bin Ranges
BB.range.Amp = linspace(-2,1.5,20);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
BB.range.segDur = linspace(2,3.75,24); %0:100:1800; % Group: 0:100:1800; single: 25:150:1800
BB.range.AmpPrc = 0:5:100;

% Setup Plot Limits
BB.plot.lims.burfreq = [0 30];
BB.plot.lims.PLV = [-100 100]; %[0 0.35];
BB.plot.lims.PLVun = [0 0.6];
BB.plot.lims.wPLV = [-10 10];
BB.plot.lims.Amp =  [-2 2];
BB.plot.lims.wAmpPrc =  [-2 6];
BB.plot.lims.Dur = [1.8 3.6];

% Compute  Burst Amplitude/Duration Statistics
BB = compute_BetaBurstAmpDurStats_v2(R,BB,[],0);

% Plot Burst Properties
subplot(1,2,1); 
[h,l] = plotBurstAmplitudeHistogram(R,BB);
delete(l);
p = get(gca,'Children');
set(gca,'Children',p([2 4 3 1])); 

subplot(1,2,2); 
plotBurstDurationHistogram(R,BB);
p = get(gca,'Children');
set(gca,'Children',p([2 4 3 1])); 

set(gcf,'Position',[631   672   784   305])

