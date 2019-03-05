function h = plotModelComparisonMatrix(R,compStruc,shortlab,type)
% TlnK = 2.*log(max(pmod)./pmod);
if type == 1 % model probs
    F = compStruc.pmod; %-log10(1-compStruc.pmod);
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    F = -log10(1-F);
    yl = 'Posterior Model Probability';
elseif type == 2 % KL divergences
    F = compStruc.KL;
    F = reshape(F,sqrt(numel(F)),sqrt(numel(F)));
    
    yl = 'KL Divergence';
elseif type == 3 % Score
    F1 = compStruc.KL;
    F1 = reshape(F1,sqrt(numel(F1)),sqrt(numel(F1)));
    F1 = log10(1./F1);
%     for i = 1:3
%         F1(i,:) = F1(i,:)./max(F1(i,:));
%     end
    F2 = compStruc.pmod;
    F2 = reshape(F2,sqrt(numel(F2)),sqrt(numel(F2)));
    F2 = -log10(1-F2);
    F = F2-F1;
    yl = 'Accuracy-Divergence';
end

h = imagesc(F)
[x,y] = meshgrid(1:size(F,1),1:size(F,2));
text(x(:),y(:),num2str(F(:),'%0.2f'),'HorizontalAlignment','center')
cmap = brewermap(128,'*RdYlBu');
colormap(cmap);
a = gca;
set(a,'YDir','rev') %'XTick',[],'YTick',[],
a.XTick = 1:size(x,1);
a.YTick = 1:size(x,1);
xlabel('Model to be Fit');
ylabel('Model Data');
title(yl);