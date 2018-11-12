function BB = BA_analysis_061118(R,BB)
%% THIS IS AN ISSUE AS OTHER SCRIPT ONLY WRITES BETA BURSTS TO STN/CORTEX
for cond = 1:size(R.condname)
    cmap = linspecer(4); %brewermap(4,'Set1');
    cmap = cmap(end:-1:1,:);
    onsetT = [-500 500];
    BEpoch = []; REpoch = [];
    for i = 1:numel(BB.segTInds{cond})
        epochdef = [BB.consecSegs{cond}{BB.consecSegsInds{cond}(i)}(1)+ floor((onsetT(1)/1e3)*BB.fsamp):BB.consecSegs{cond}{BB.consecSegsInds{cond}(i)}(1) + floor((onsetT(2)/1e3)*BB.fsamp)];
        if epochdef(end)<size(BB.A{1},2) && epochdef(1) > 0
            BEpoch(:,:,i) = BB.A{cond}(:,epochdef);
            REpoch(:,:,i) = BB.raw{cond}(:,epochdef);
        end
    end
    TEpoch = linspace(onsetT(1),onsetT(2),size(epochdef,2));
    meanEnv = squeeze(mean(BEpoch,3)); stdEnv = squeeze(std(BEpoch,[],3)); %./sqrt(size(BEpoch,3)));
    meanRaw = squeeze(mean(REpoch,3));
    [x ind] = max(meanEnv,[],2)
    
    
    subplot(1,3,cond)
    for i = 1:size(BEpoch,1)
        XY = 1*meanEnv(i,:);
        XYp = XY-mean(XY);
        XY = 10 + ((XYp)-(i*2));
        lhan((2*i)-1) = plot(TEpoch,XY','color',cmap(i,:),'LineWidth',2); hold on
        [lp hp] = boundedline(TEpoch,XY',stdEnv(i,:)')
        lp.Color = cmap(i,:);
        hp.FaceColor = cmap(i,:); hp.FaceAlpha = 0.5;
%         XY = 1.2*(10.^XYp).*meanRaw(i,:);
            XY = 3.5*meanRaw(i,:);
        XY = XY - mean(XY);
        XY = 10 + ((XY)-(i*2));
        plot(TEpoch,XY','color',cmap(i,:),'LineWidth',1.5)
        lhan(2*i) = plot(repmat(TEpoch(ind(i)),2,1),[1.5 500],'LineStyle','--','color',cmap(i,:),'LineWidth',1);
        splay = TEpoch(ind(i))+(10*(TEpoch(ind(i))-median(TEpoch(ind))));
        plot([splay TEpoch(ind(i))],[0.5 1.5],'LineStyle','--','color',cmap(i,:),'LineWidth',1)
        text(splay,0.4,sprintf('%0.1f ms',TEpoch(ind(i))),'color',cmap(i,:))
        
        legn{(2*i)-1} = R.chsim_name{i};
        legn{2*i} = [R.chsim_name{i} ' max'];
    end
    % legend(lhan(1:2:8),legn(1:2:8))
    
    xlabel('STN Beta Onset Time (ms)'); ylabel('Average Amplitude (a.u.)'); ylim([0 10]); xlim(onsetT)
    set(gca,'YTickLabel',[]);
    title('Sequential Beta Onset in Model 5')
end
set(gcf,'Position',[680          85         419        1013])