function parSweepPlot(R,parsweep)
colormap jet
set(gcf,'color','w');
for i = 1:numel(R.chsim_name)
subplot(2,3,1)
imagesc2(parsweep.Rlist,parsweep.Qlist,squeeze(parsweep.betaPowBank(2,:,:)))
set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
xlabel('STN->GPe','FontSize',12)
ylabel('GPe-|STN','FontSize',12)

axis square
title('Membrane surface area at equilibrium','FontSize',12,'fontweight','bold')
c = colorbar;
ylabel(c,'log Membrane surface area(cm^2)','FontSize',12)
end