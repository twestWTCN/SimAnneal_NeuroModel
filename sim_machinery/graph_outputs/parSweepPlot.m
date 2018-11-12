function parSweepPlot(R,parsweep,cmap)
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery\graph_outputs\cmap_ABC.mat')
cmap = linspecer(32);
colormap(cmap)
set(gcf,'color','w');

Xr = min(parsweep.Rlist):.001:max(parsweep.Rlist);
Xq = min(parsweep.Qlist):.001:max(parsweep.Qlist);
[X,Y] = meshgrid(Xr,Xq);
% V = interp2(parsweep.Rlist,parsweep.Qlist,squeeze(parsweep.betaPowBank(i,:,:)),X,Y,'spline');
V = interp2(parsweep.Rlist,parsweep.Qlist,parsweep.plot.feat,X,Y,'spline');
imagesc2(Xr,Xq,V)
set(gca,'XAxisLocation','bottom')
set(gca,'YDir','normal')
