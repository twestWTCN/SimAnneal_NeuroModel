function parSweepPlot(R,parsweep)
colormap jet
set(gcf,'color','w');
for i = 1:numel(R.chsim_name)
subplot(2,3,i)
Xr = min(parsweep.Rlist):.01:max(parsweep.Rlist);
Xq = min(parsweep.Qlist):.01:max(parsweep.Qlist);
[X,Y] = meshgrid(Xr,Xq);
V = interp2(parsweep.Rlist,parsweep.Qlist,squeeze(parsweep.betaPowBank(i,:,:)),X,Y,'cubic');
imagesc2(Xr,Xq,V)
set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
xlabel('STN->GPe','FontSize',10)
ylabel('GPe-|STN','FontSize',10)
axis square
title([R.chsim_name{i}],'FontSize',10,'fontweight','bold')
% c = colorbar;
end
set(gcf,'Position',[610.0000  157.5000  950.0000  635.5000])
annotation(gcf,'textbox',...
    [0.152052631578947 0.962399596025029 0.745842105263159 0.0299837925445702],...
    'String',{'Parameter Sweep of GPe/STN Coupling and Resultant Local Beta Power'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontWeight','bold',...
    'FitBoxToText','off');
