function [hl,hp,hs] =  plotBoundedScatter(pirange,Yb,Ybz,Y,X,cmap)
    [hl, hp] = boundedline(pirange,Yb,Ybz);
     hl.Color = cmap; hp.FaceColor = cmap; hp.FaceAlpha = 0.3; hl.LineWidth = 2;
    hold on
    hs = scatter(Y,X,45);
    hs.MarkerFaceColor = cmap;
    hs.MarkerEdgeColor = 'none';
    grid on
xlim([-pi pi])