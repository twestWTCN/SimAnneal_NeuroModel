function plotTLPhaseLocking(TL,cond,feat,ip)
% phslipN = sum(phslip,1);
% [B,I] = sort(BB.segDur{1});
% phslip = phslip(:,I)
%  subplot(1,2,1)
% imagesc(TL.epochT,1:size(phslip,2),phslip')
% phslipProb =sum(phslip,2)/size(phslip,2);
%  xlim([-300 300])

for L = [1] % M2 STR GPe
    %     if L == 1;
    %             title(R.condname{cond})
    %     end
    
    %     subplot(size(TL.(feat){cond},1),3,plist(ip,L))
    subplot(1,3,ip)
    for dur = 1:2
    phi = squeeze(TL.phi{cond}([L 4],:,TL.durFlag{cond}(:,dur)));
    
    winsize = floor(size(phi,2)/15);
    swInds = slideWindow(1:size(phi,2),winsize,floor(0.85*winsize));
    swInds(swInds==0) = NaN;
    swTEpoch = []; psProb = [];
    for p = 1:size(swInds,2)-1
        swTEpoch(p) = median(TL.epochT(swInds(:,p)));
        dphi = squeeze(phi(1,swInds(:,p),:)) - squeeze(phi(2,swInds(:,p),:));
        PLV(:,p) =  abs(mean(exp(-1i*dphi),2));
    end
    
    [hl, hp] = boundedline(swTEpoch,median(PLV,1),std(PLV)); %-(L.*shift)
    hp.FaceColor = TL.struccmap(L,:).*0.75;
    hp.FaceAlpha = 0.75;
    hl.Color = TL.struccmap(L,:).*0.75;
    hl.LineWidth = 1;
    if dur == 2
    hl.Color = TL.struccmap(L,:);
    hp.FaceColor = TL.struccmap(L,:);
    hl.LineStyle = '--';
    hp.LineStyle = '--';
    end
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
end
