close all
CON = 1;
%     load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'])
%

uq = prctile(BB.segDur{cond},75);
lq = prctile(BB.segDur{cond},25);
twin = linspace(-250/2000,250/2000,size(BB.segEnv{6},2))*1000
cmap = linspecer(4);
for i = 1:4
    %     plot(twin,mean(BB.segEnv{6}(i,:,BB.segDur{cond}<lq),3),'color',cmap(i,:))
    boundedline(twin,squeeze(mean(BB.segEnv{6}(i,:,BB.segDur{cond}>uq),3)),squeeze(std(BB.segEnv{6}(i,:,BB.segDur{cond}>uq),[],3)./sqrt(size(BB.segEnv{6},3))),'cmap',cmap(i,:));
    hold on
    %     a(i) =  plot(twin,mean(BB.segEnv{6}(i,:,BB.segDur{cond}>uq),3),'color',cmap(i,:),'LineWidth',2)
    [a(i) b] = boundedline(twin,squeeze(mean(BB.segEnv{6}(i,:,BB.segDur{cond}<lq),3)),squeeze(std(BB.segEnv{6}(i,:,BB.segDur{cond}>uq),[],3)./sqrt(size(BB.segEnv{6},3))),'cmap',cmap(i,:));
    a(i).LineWidth = 2;
end

legend(a,R.chloc_name)

figure
for  i = 1:4
    
    acf = [];
    for j = find(BB.segDur{cond}<lq)
        [acf(:,j),lags] = autocorr(squeeze(BB.segEnv{6}(i,:,j)),'NumLags',500);
    end
    plot(lags,mean(acf,2),'color',cmap(i,:))
    
    acf = [];
    for j = find(BB.segDur{cond}>uq)
        [acf(:,j),lags,bounds] = autocorr(squeeze(BB.segEnv{6}(i,:,j)),'NumLags',500);
    end
    hold on
    plot(lags,mean(acf,2),'color',cmap(i,:),'LineWidth',2)
end


figure
chcmb = nchoosek(1:4, 2);
cmap = linspecer(6)
for  i = size(chcmb,1)
    
    acf = [];
    for j = find(BB.segDur{cond}<lq)
        [acf(:,j),lags] = xcorr(squeeze(BB.segEnv{6}(chcmb(i,1),:,j)),squeeze(BB.segEnv{6}(chcmb(i,2),:,j)),'coeff');
    end
    plot(lags,mean(acf,2),'color',cmap(i,:))
    
    acf = [];
    for j = find(BB.segDur{cond}>uq)
        [acf(:,j),lags] = xcorr(squeeze(BB.segEnv{6}(chcmb(i,1),:,j)),squeeze(BB.segEnv{6}(chcmb(i,2),:,j)),'coeff');
    end
    hold on
    plot(lags,mean(acf,2),'color',cmap(i,:),'LineWidth',2)
end



X =  squeeze(mean(BB.segEnv{6}(4,1:250,:),2))
Y = BB.segDur{6};
