function plotModelComparisonMatrix(R,compStruc,shortlab,type)
% TlnK = 2.*log(max(pmod)./pmod);
if type == 1 % model probs
    F = -log10(1-compStruc.pmod);
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    yl = '-log_{10} P(M|D)';
elseif type == 2 % KL divergences
    F = compStruc.KL;
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    
    yl = 'KL Divergence';
elseif type == 3
    F1 = compStruc.KL;
    F1 = reshape(F1,sqrt(numel(F1)),sqrt(numel(F1)));
    for i = 1:3
        F1(i,:) = F1(i,:)./max(F1(i,:));
    end
    F2 = compStruc.pmod;
    F2 = reshape(F2,sqrt(numel(F2)),sqrt(numel(F2)));
    F = F2-F1;
    yl = 'Score';
end

imagesc(F)
[x,y] = meshgrid(1:size(F,1),1:size(F,2));
text(x(:),y(:),num2str(F(:)),'HorizontalAlignment','center')
colormap(R.plot.cmap);
set(gca,'YDir','normal') %'XTick',[],'YTick',[],

a.XTickLabel = shortlab;
xlabel('Sim Model');
ylabel('Fit Model');
title(yl);