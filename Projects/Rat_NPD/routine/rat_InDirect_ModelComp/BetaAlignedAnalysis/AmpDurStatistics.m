function BB = AmpDurStatistics(R,BB)
close all
for CON = 1:4
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'])
    figure
    % Setup Bin Ranges
    BB.range.Amp = linspace(0,20,20);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
    BB.range.Amplr = linspace(0,20,10);% 1:0.25:5; %0:5:120; % Group: 0:3:120; single: 0:5:80
    % BB.range.Dur = linspace(50,1400,20);
    BB.range.Dur = linspace(log10(20),log10(1800),20);
    BB.range.Durlr = linspace(log10(20),log10(1800),10);
    
    BB.range.segDur = linspace(1.5,3.2,24); %0:100:1800; % Group: 0:100:1800; single: 25:150:1800
    BB.range.AmpPrc = 0:5:100;
    
    % Setup Plot Limits
    BB.plot.lims.burfreq = [0 30; 0 30];
    BB.plot.lims.PLV = [-100 100]; %[0 0.35];
    BB.plot.lims.PLVun = [0 0.6];
    BB.plot.lims.wPLV = [-10 10];
    BB.plot.lims.Amp =  [0 20];
    BB.plot.lims.wAmpPrc =  [-2 6];
    BB.plot.lims.Dur = log10([15 1800]);
    BB.plot.durlogflag = 1;
    % Compute  Burst Amplitude/Duration Statistics
    % for i = 1:5; F(i) = figure; end
    
    BB = computeBetaBurstAmpDurStats_v2(R,BB);
    R.condcmap = linspecer(30);
%     R.condname([1 6 8 15 17 20]) =  {'1% M2->STN','Fitted','150% M2->STN','1% STN ->GPe','Fitted','150% STN->GPe'};
    R.condnames =  R.condname;
    
    % Plot Burst Properties
    
    %     if CON == 1
    condsel = [20 18 16]; % conditions to be selected
    R.condcmap(condsel,:) =  linspecer(3);
    ip(:,1) = [1 2 3]; % panel indice for plotting
    %     else
    %         condsel = [15 17 20];
    %         ip(:,2) = [2 4 6];
    %     end
    %     BB = compute_BurstThreshold(R,BB,condsel,0);
    %     R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STR->GPe','150% STR->GPe'};
    
    subplot(2,3,1);
    [h,l] = plotBurstAmplitudeHistogram(R,BB,condsel);
    ylim(BB.plot.lims.burfreq(1,:))
    p = get(gca,'Children');
    %     if CON == 1
    %         set(gca,'Children',p([3 2 1]));
    %     else
    %         set(gca,'Children',p([3 2 1]));
    %     end
    subplot(2,3,2);
    plotBurstDurationHistogram(R,BB,condsel);
    ylim(BB.plot.lims.burfreq(2,:))
    p = get(gca,'Children');
    p(1).FaceAlpha = 0.85;
    %     if CON == 1
    %         set(gca,'Children',p([3 2 1]));
    %     else
    %         set(gca,'Children',p([3 2 1]));
    %     end
    subplot(2,3,3);
    plotBurstAmpDurScatter(R,BB,condsel)
    %     p = get(gca,'Children');
    %     set(gca,'Children',p([1 3 5 2 4 6]));
    
    
    %     cmap = brewermap(30,'Reds');
    %     cntlevs = 0.5:0.5:15;
    %     [h,l] = plotBurstAmplitudeDurationHistogram(R,BB,condsel,cmap,cntlevs);
    %     cb = colorbar('southoutside');
    %     cb.Position = [0.1548    0.0360    0.7589    0.0112];
    %     cb.Label.String = 'Burst Rate (min^-1)';
    %
    
    
    
    condsel = find(BB.condstat.ssAmp(1,:)<1e3);
    featlist = {'ssAmp','ssDur','ssPLV'};
    featname = {'Burst Amplitude','Burst Duration','Burst PLV'};
    yl_list = {[0 10],[0 400],[0.5 1]}
    for feat = 1:3
        subplot(2,3,feat+3)
        boundedline(log10(BB.condlist(condsel)),BB.condstat.(featlist{feat})(1,condsel),BB.condstat.(featlist{feat})(2,condsel)./2)
        
        xlabel('Connection Strength')
        ylabel(featname{feat})
        xlim([-1 0.75])
        ylim(yl_list{feat})
        grid on
    end
    
    
    
    
    set(gcf,'Position',[313    45   1300    675])
    %     figure
    %     cmap = brewermap(30,'Blues');
    %     cntlevs = 0.5:0.5:10;
    %     [h,l] = plotBurstAmplitudeDurationHistogram(R,BB,[15 17 20],cmap,cntlevs);
    %     cb = colorbar('southoutside');
    %     cb.Position = [0.1548    0.0360    0.7589    0.0112];
    %     cb.Label.String = 'Burst Rate (min^-1)';
    %     set(gcf,'Position',[680    29   336   949])
end
