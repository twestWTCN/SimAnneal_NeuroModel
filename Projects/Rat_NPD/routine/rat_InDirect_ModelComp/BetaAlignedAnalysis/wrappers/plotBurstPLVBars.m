close all

condSel = [1:20];
stats_PPC = [];
p = 0;
splitlab = {'Duration < 25%','Duration > 75%'};
for cond = condSel
    IndsByDur{1} = BB.segDur{cond}<prctile(BB.segDur{cond},25);
    IndsByDur{2} = BB.segDur{cond}>prctile(BB.segDur{cond},75);
    p = p + 1;
    
    stats_PPC(:,p,1) = [median(BB.segPPC{cond}(IndsByDur{1})) median(BB.segPPC{cond}(IndsByDur{2}))];
    stats_PPC(:,p,2) = [iqr(BB.segPPC{cond}(IndsByDur{1}))./sqrt(sum(IndsByDur{1}))
        iqr(BB.segPPC{cond}(IndsByDur{2}))./sqrt(sum(IndsByDur{2}))];
    
    stats_amp(:,p,1) = [median(BB.segAmpPrc{cond}(IndsByDur{1})) median(BB.segAmpPrc{cond}(IndsByDur{2}))];
    stats_amp(:,p,2) = [iqr(BB.segAmpPrc{cond}(IndsByDur{1}))./sqrt(sum(IndsByDur{1}))
        iqr(BB.segAmpPrc{cond}(IndsByDur{2}))./sqrt(sum(IndsByDur{2}))];
    
    
    % BB.range.RP = linspace(-pi,pi,8);
    % [BB.Amp_binRP(:,:,cond),BB.Amp_binRP_data{cond}] = binDatabyRange(BB.segRP{cond}(IndsByDur{1}),BB.range.RP,BB.segAmp{cond}(IndsByDur{1}));
    % [BB.Amp_binRP(:,:,cond),BB.Amp_binRP_data{cond}] = binDatabyRange(BB.segRP{cond}(IndsByDur{1}),BB.range.RP,BB.segAmp{cond}(IndsByDur{1}));
    
end
subplot(2,1,1)
H = bar(1:20,stats_PPC(:,:,1)');
hold on

a= gca;
a.XTickLabel = R.condname(condSel);
a.YLabel.String = 'M2/STR Pairwise Phase Consistency';
es = stats_PPC(:,:,2)';
em = stats_PPC(:,:,1)';

% sp = [0.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15 4.85 5.15];
% E = errorbar(sp,reshape(em',1,10),reshape(es',1,10));
% E.LineStyle = 'none';

subplot(2,1,2)
H = bar(1:20,stats_amp(:,:,1)');
hold on

a= gca;
a.XTickLabel = R.condname(condSel);
a.YLabel.String = 'Median Burst Amplitude';
es = stats_amp(:,:,2)';
em = stats_amp(:,:,1)';

% sp = [0.85 1.15 1.85 2.15 2.85 3.15 3.85 4.15 4.85 5.15];
% E = errorbar(sp,reshape(em',1,10),reshape(es',1,10));
% E.LineStyle = 'none';
% legend(H,splitlab)

figure
BB.range.RPmid = binEdge2Mid(BB.range.RP);
condSelG = [1:10; 11:20];
for i = 1:2
    if i == 1
    R.condcmap = brewermap(10,'Reds');
    else
    R.condcmap(11:20,:) = brewermap(10,'Blues');
    end
    subplot(2,2,i)
    for cond = condSelG(i,:)
        plot(BB.range.RPmid,squeeze(BB.Dur_binRP(1:end-1,1,cond)),'color',R.condcmap(cond,:),'LineWidth',2)
        hold on
    end
    xlim([-pi pi])
    ylabel('Median Burst Duration (ms)')
    xlabel('Relative Phase (\phi_{M2}-\phi_{STN})')
    
    subplot(2,2,2+i)
    for cond = condSelG(i,:)
        plot(BB.range.RPmid,squeeze(BB.AmpPrc_binRP(1:end-1,1,cond)),'color',R.condcmap(cond,:),'LineWidth',2)
        hold on
    end
    xlim([-pi pi])
    ylabel('Median Burst Amplitude (ms)')
    xlabel('Relative Phase (\phi_{M2}-\phi_{STN})')    
end