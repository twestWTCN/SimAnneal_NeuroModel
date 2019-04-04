function [hl,hp,hs] =  plotBoundedScatter(Xrange,Yb,Ybz,Y,X,cmap)
    [hl, hp] = boundedline(Xrange,Yb,Ybz);
     hl.Color = cmap; hp.FaceColor = cmap; hp.FaceAlpha = 0.3; hl.LineWidth = 2;
    hold on
    hs = scatter(Y,X,35);
    hs.MarkerFaceColor = cmap;
    hs.MarkerEdgeColor = 'none';
    grid on
