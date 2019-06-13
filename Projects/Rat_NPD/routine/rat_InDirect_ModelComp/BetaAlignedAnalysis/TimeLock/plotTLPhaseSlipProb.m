function plotTLPhaseSlipProb(TL,cond,feat,ip)
% phslipN = sum(phslip,1);
% [B,I] = sort(BB.segDur{1});
% phslip = phslip(:,I)
%  subplot(1,2,1)
% imagesc(TL.epochT,1:size(phslip,2),phslip')
% phslipProb =sum(phslip,2)/size(phslip,2);
%  xlim([-300 300])
plist(1,:) = [1 4 7 10];
plist(2,:) = [2 5 8 11];
plist(3,:) = [3 6 9 12];

for L = 1:size(TL.(feat){cond},1)
%     if L == 1;
%             title(R.condname{cond})
%     end
          
%     subplot(size(TL.(feat){cond},1),3,plist(ip,L))
    subplot(1,3,ip)
    
    PS = squeeze(TL.(feat){cond}(L,:,:));
    phslip1 = double(PS>prctile(PS(:),75));
    
    swInds = slideWindow(1:size(phslip1,1), floor(size(phslip1,1)/300));
    swInds(swInds==0) = NaN;
    swTEpoch = []; psProb = [];
    for i = 1:size(swInds,2)-1
            swTEpoch(i) = median(TL.epochT(swInds(:,i)));
            psProb(i,1) = sum(sum(phslip1(swInds(:,i),:),2))/numel(phslip1(swInds(:,i),:));
            psProb(i,2) = sqrt(sum(sum(phslip1(swInds(:,i),:),2)))/numel(phslip1(swInds(:,i),:));
%             psProb(i,3) = sum(sum(phslip3(swInds(:,i),:),2))/numel(phslip3(swInds(:,i),:));
    end
    psProb(1,:) = nan(1,2);
%     psProb(:,1) = sum(phslip,2)/size(phslip,2);
%     
%     phslip = double(squeeze(TL.(feat){cond}(L,:,:))>eps*0.99);
%     psProb(:,2) = sum(phslip,2)/size(phslip,2);
%     phslip = double(squeeze(TL.(feat){cond}(L,:,:))>eps*1.01);
%     psProb(:,3) = sum(phslip,2)/size(phslip,2);
    
    [hl, hp] = boundedline(swTEpoch,psProb(:,1),psProb(:,2)); %-(L.*shift)
    hp.FaceColor = TL.struccmap(L,:);
    hp.FaceAlpha = 0.75;
    hl.Color = TL.struccmap(L,:);
    hl.LineWidth = 1;
    %     p(L) = plot(TL.epochT,mean(TL.(feat){cond}(L,:,:),3)-(L.*shift),'color',TL.struccmap(L,:));
    %     p(L).LineWidth = 2;
    hold on
%     set(gca, 'YScale', 'log')
%     if L>2
%     ylim([5e-3 0.5])
%     elseif L == 1
%     ylim([0.05 0.3])
%     elseif L == 2
%         ylim([0.1 0.75])
%     end

end
 