function plotTMS_MEPStats(feat,cmap,lwid)
% xrange = linspace(-pi,pi,8);
% feat_emp{1} = [binEdge2Mid(xrange);binStatsRS_phiAmp'];
% feat_emp{2} = [binEdge2Mid(xrange);binStatsRS_phiAmpCoV'];
% feat_emp{3} = [binEdge2Mid(xrange);binStatsRS_phiLag'];
% feat_emp{4} = [binEdge2Mid(xrange);binStatsRS_phiLagCoV'];
% xrange = 0.5:0.75:35;
% feat_emp{5} = [binEdge2Mid(xrange);binStats_ampLag(:,3)'];
% feat_emp{6} = [lagX;lagPDF];
% feat_emp{7} = [TMS_phase; MEP_amp; MEP_onset];


% TMS_phase = feat{end}(1,:);
% MEP_amp = feat{end}(2,:);
% MEP_onset = feat{end}(3,:);

subplot(4,1,1)
% a(1) = scatter(TMS_phase(1:10:end),MEP_onset(1:10:end),15);
% hold on
% a(1).CData = cmap;
% yyaxis right
% a(2) =scatter(TMS_phase(1:10:end),MEP_amp(1:10:end),15);
% a(2).CData = cmap;
% a(2).MarkerFaceColor = cmap;
% a(2).MarkerFaceAlpha = 0.5;
yyaxis left
p = plot(feat{1}(1,:),feat{1}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','-');
hold on
yyaxis right
p = plot(feat{3}(1,:),feat{3}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','--'); 
grid on
xlabel('EEG Phase at TMS onset')
yyaxis left
ylabel('MEP Amplitude');
yyaxis right
ylabel('MEP Latency')

subplot(4,1,2)
yyaxis left
p = plot(feat{2}(1,:),feat{2}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','-');
hold on
% yyaxis right
% p = plot(feat{4}(1,:),feat{4}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','--'); 

grid on
xlabel('EEG Phase at TMS onset')
yyaxis left
ylabel('MEP Amplitude CoV ');
% yyaxis right
% ylabel('MEP Latency CoV')

subplot(4,1,3)
% a(3) = scatter(MEP_amp(1:10:end),MEP_onset(1:10:end));
% a(3).CData = cmap;
% a(3).MarkerFaceColor = cmap;
% a(3).MarkerFaceAlpha = 0.58;
p = plot(feat{5}(1,:),feat{5}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','--'); 
hold on
xlabel('MEP Z Amplitude')
ylabel('MEP Mean CorrectedAmplitude');

subplot(4,1,4)
% h = histogram(MEP_onset,12,'Normalization','pdf');
% h.FaceColor = cmap;
% h.FaceAlpha = 0.8;
plot(feat{6}(1,:),feat{6}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','--')
hold on
plot(feat{7}(1,:),feat{7}(2,:),'LineWidth',lwid,'color',cmap,'LineStyle','-')
ylabel('P(x)');
xlabel('MEP Latency')
