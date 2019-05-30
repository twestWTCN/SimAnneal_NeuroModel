function p = plotTLTimeEvolutions(TL,cond,feat,shift)
if nargin<4
shift = TL.(feat){cond};
shift = std(shift(:));
end
for L = 1:size(TL.(feat){cond},1)
    p(L) = plot(TL.epochT,mean(TL.(feat){cond}(L,:,:),3)-(L.*shift),'color',TL.struccmap(L,:));
    p(L).LineWidth = 2;
    hold on
end
