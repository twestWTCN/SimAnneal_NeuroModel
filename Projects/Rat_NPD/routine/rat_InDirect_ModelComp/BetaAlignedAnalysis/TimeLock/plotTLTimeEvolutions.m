function p = plotTLTimeEvolutions(TL,cond,feat,shift)
if nargin<4
    shift = TL.(feat){cond};
    shift = std(shift(:));
end
for L = 1:size(TL.(feat){cond},1)
    
    [hl, hp] = boundedline(TL.epochT,mean(TL.(feat){cond}(L,:,:),3)-(L.*shift),2.326.*std(TL.(feat){cond}(L,:,:),[],3)./sqrt(size(TL.(feat){cond},3)));
    hp.FaceColor = TL.struccmap(L,:);
    hp.FaceAlpha = 0.75;
    hl.Color = TL.struccmap(L,:);
    hl.LineWidth = 1;
    %     p(L) = plot(TL.epochT,mean(TL.(feat){cond}(L,:,:),3)-(L.*shift),'color',TL.struccmap(L,:));
    %     p(L).LineWidth = 2;
    hold on
end
 