function plotSimulatedPRC(PRC,condsel)

for i = condsel
    subplot(1,2,2)
    a = scatter(PRC.impPhi_cS{i},PRC.impdPhi_cS{i});
    a.MarkerEdgeColor = PRC.condcmap(i,:);
    a.MarkerFaceColor = PRC.condcmap(i,:);
    a.MarkerFaceAlpha = 0.5;
    a.LineWidth = 1;
    hold on
    
    [xq yq R2] = sinfit(PRC.impPhi_cS{i},PRC.impdPhi_cS{i},50,0);
    b = plot(xq,yq);
    b.Color = PRC.condcmap(i,:).*0.8;
    b.LineWidth = 2;
    xlabel('Impulse Phase \phi');
    ylabel('Change in phase \Delta\phi')
    
    subplot(1,2,1)
    a = scatter(PRC.impPhi_cS{i},PRC.impdA_cS{i});
    a.MarkerEdgeColor = PRC.condcmap(i,:);
    a.MarkerFaceColor = PRC.condcmap(i,:);
    a.MarkerFaceAlpha = 0.5;
    a.LineWidth = 1;
    hold on
    
    [xq yq R2] = sinfit(PRC.impPhi_cS{i},PRC.impdA_cS{i},25,0);
    b = plot(xq,yq);
    b.Color = PRC.condcmap(i,:).*0.8;
    b.LineWidth = 2.5;
    xlabel('Impulse Phase \phi');
    ylabel('Change in Amplitude')
end

subplot(1,2,1)
grid on
p = get(gca,'Children');
set(gca,'Children',p([1 3 5 2 4 6]))
hl = findobj(p,'type','line');
title('Phase Response Curve- Amplitude')

subplot(1,2,2)
grid on
p = get(gca,'Children');
set(gca,'Children',p([1 3 5 2 4 6]))
ylim([-1 1])
hl = findobj(p,'type','line');
title('Phase Response Curve- Phase Shift')

legend(hl,PRC.condname(condsel),'Location','SouthEast')