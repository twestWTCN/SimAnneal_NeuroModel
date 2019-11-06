function plotSimulatedPRC(PRC,condsel,conStren)
cmap = brewermap(18,'Spectral');
for i = condsel
    B_sel = PRC.impBetaLev{i}>0;
    subplot(1,3,2)
    a = scatter(PRC.impPhi_cS{i}(B_sel),PRC.impdPhi_cS{i}(B_sel),PRC.impBetaLev{i}(B_sel)-min(PRC.impBetaLev{i}(B_sel))+1);
    a.MarkerEdgeColor = PRC.condcmap(i,:);
    a.MarkerFaceColor = PRC.condcmap(i,:);
    a.MarkerFaceAlpha = 0.5;
    a.LineWidth = 1;
    hold on
    
    [xq yq R2] = sinfit(PRC.impPhi_cS{i}(B_sel),PRC.impdPhi_cS{i}(B_sel),50,0);
    b = plot(xq,yq);
    b.Color = PRC.condcmap(i,:).*0.8;
    b.LineWidth = 2;
    xlabel('Impulse Phase \phi');
    ylabel('Change in phase \Delta\phi')
    
    subplot(1,3,1)
    a = scatter(PRC.impPhi_cS{i}(B_sel),PRC.impdA_cS{i}(B_sel));
    a.MarkerEdgeColor = PRC.condcmap(i,:);
    a.MarkerFaceColor = PRC.condcmap(i,:);
    a.MarkerFaceAlpha = 0.5;
    a.LineWidth = 1;
    hold on
    
    [xq yq R2] = sinfit(PRC.impPhi_cS{i}(B_sel),PRC.impdA_cS{i}(B_sel),25,0);
    b = plot(xq,yq);
    b.Color = PRC.condcmap(i,:).*0.8;
    b.LineWidth = 2.5;
    xlabel('Impulse Phase \phi');
    ylabel('Change in Amplitude')
    
    
    subplot(1,3,3)
   binRP = binDatabyRange(PRC.impPhi_cS{i}(B_sel),-pi:pi/6:pi,PRC.impdPhi_cS{i}(B_sel),'number');
    a = scatter(binEdge2Mid(-pi:pi/6:pi),binRP(1:end-1,4));
    a.MarkerEdgeColor = PRC.condcmap(i,:);
    a.MarkerFaceColor = PRC.condcmap(i,:);
    a.MarkerFaceAlpha = 0.5;
    hold on
    
    a = plot(binEdge2Mid(-pi:pi/6:pi),binRP(1:end-1,4));
    a.Color = PRC.condcmap(i,:);
    a.LineWidth = 2;
    
    
end


%% Summary STats of PRC
figure(3)

subplot(2,1,1)
% For Amplitudes
for i = 1:size(PRC.impBetaLev,2)
    B_sel = PRC.impBetaLev{i}>0;
    mindPhi(i) = min(PRC.impdA_cS{i}(B_sel));
    maxdPhi(i) = max(PRC.impdA_cS{i}(B_sel));
    randPhi(i) = range(PRC.impdA_cS{i}(B_sel));
end
[p sc] = plotPRCSumStats(conStren,maxdPhi,mindPhi,randPhi,condsel,condcmap);

set(p,'Color',cmap(end,:),'LineWidth',2)
ylabel('% Change in STN Beta Amplitude')
xlabel('Connection Strength (% of fitted)')
grid on
legend(p,'PRC Max','PRC Min','PRC Range')
xlim([5 250])


subplot(2,1,2)
% For Angles
for i = 1:size(PRC.impBetaLev,2)
    B_sel = PRC.impBetaLev{i}>0;
    mindPhi(i) = min(PRC.impdPhi_cS{i}(B_sel));
    maxdPhi(i) = max(PRC.impdPhi_cS{i}(B_sel));
    randPhi(i) = range(PRC.impdPhi_cS{i}(B_sel));
end
[p sc] = plotPRCSumStats(conStren,maxdPhi,mindPhi,randPhi,condsel,condcmap);




% subplot(1,2,1)
% grid on
% p = get(gca,'Children');
% set(gca,'Children',p([1 3 5 2 4 6]))
% hl = findobj(p,'type','line');
% title('Phase Response Curve- Amplitude')
% 
% subplot(1,2,2)
% grid on
% p = get(gca,'Children');
% set(gca,'Children',p([1 3 5 2 4 6]))
% ylim([-1 1])
% hl = findobj(p,'type','line');
% title('Phase Response Curve- Phase Shift')
% 
% legend(hl,PRC.condname(condsel),'Location','SouthEast')