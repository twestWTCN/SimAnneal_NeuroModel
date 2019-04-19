function BB = AmpDurStatistics(R,BB)
% Setup Bin Ranges
BB.range.Amp = linspace(-2,-0.5,20);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
BB.range.segDur = linspace(1.5,3.2,24); %0:100:1800; % Group: 0:100:1800; single: 25:150:1800
BB.range.AmpPrc = 0:5:100;

% Setup Plot Limits
BB.plot.lims.burfreq = [0 6; 0 15];
BB.plot.lims.PLV = [-100 100]; %[0 0.35];
BB.plot.lims.PLVun = [0 0.6];
BB.plot.lims.wPLV = [-10 10];
BB.plot.lims.Amp =  [-2 -0.5];
BB.plot.lims.wAmpPrc =  [-2 6];
BB.plot.lims.Dur = [1.8 3.2];

% Compute  Burst Amplitude/Duration Statistics
% for i = 1:5; F(i) = figure; end
BB = computeBetaBurstAmpDurStats_v2(R,BB);

% Plot Burst Properties
for i = 1:2
    if i == 1
        condsel = [1 2 3]; % conditions to be selected
        ip(:,1) = [1 3 5]; % panel indice for plotting
    else
        condsel = [1 4 5];
        ip(:,2) = [2 4 6];
    end
    
    BB = compute_BurstThreshold(R,BB,condsel,0);
    
    subplot(3,2,ip(1,i));
    [h,l] = plotBurstAmplitudeHistogram(R,BB,condsel);
    ylim(BB.plot.lims.burfreq(i,:))
    p = get(gca,'Children');
    set(gca,'Children',p([2 3 1]));
    
    subplot(3,2,ip(2,i));
    plotBurstDurationHistogram(R,BB,condsel);
    ylim(BB.plot.lims.burfreq(i,:))
    p = get(gca,'Children');
    set(gca,'Children',p([2 3 1]));
    
    subplot(3,2,ip(3,i));
    plotBurstAmpDurScatter(R,BB,condsel)    
    p = get(gca,'Children');
    set(gca,'Children',p([1 3 5 2 4 6]));
end
set(gcf,'Position',[680 104 1024 800])

