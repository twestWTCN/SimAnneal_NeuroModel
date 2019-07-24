function plotSweepSpectraWrapper(R)
close all
% load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat_F1.mat'],'feat_HD','feat_STR_GPe')
cmap = brewermap(30,'Spectral');
R.CONnames = {'M2 -> STN','STR -| GPe','GPe -| STN','STN -> GPe'};
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(30,'Reds');
% cmap1(22,:) = [0 0 0];
cmap2 = brewermap(30,'Blues');
% cmap2(28,:) = [0 0 0];

figure(1)
for CON = 1:4
    load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat_F1.mat'])
    subplot(2,4,CON)
    plotSweepSpectra(R.frqz,feat,feat{6},cmap1,{R.condname{[2 1 3]}},[1 15 30],1:2:30)
    title(R.CONnames{CON})
    ylim([1e-16 1e-11])
end
% set(gcf,'Position',[684         501        1024         366])

ck_1 = logspace(-1,1,30);

% cmap = brewermap(30,'Reds');

figure(1)

for CON = 1:4
    load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat_F1.mat'])
   bpow = []; fpow = [];
    for ck = 1:numel(feat)
        [bpowr(ck) b] = max(feat{ck}(1,4,4,3,:));
        fpow(ck) =R.frqz(b);
        
    end
    
    % Scale bpow to 0
    [a zind] = min(abs(ck_1-1));
    bpow = 100*(bpowr)/bpowr(zind)
    % Remove non-physiological
    powInds = find(bpowr>1e-8);
    fpow(powInds) = nan(1,numel(powInds));
    bpow(powInds) = nan(1,numel(powInds));
    
    
    subplot(2,4,CON+4)
    br = plot(log10(ck_1(1,:)),log10(bpow),'k-');
    hold on
    Sr = scatter(log10(ck_1(1,:)),log10(bpow),50,cmap1,'filled');
    ylabel('log % of Fitted Power')
    grid on
    
    yyaxis right
    br = plot(log10(ck_1(1,:)),(fpow),':');
    hold on
    Sr = scatter(log10(ck_1(1,:)),(fpow),50,cmap2,'filled');
    Sr.Marker = 'diamond';
    ylabel('Peak Frequency (Hz)')
    xlabel('log_{10} % Connection Strength')
    title(R.CONnames{CON})
    xlim([-1 0.7]); ylim([12 25])

    
    yyaxis left
    ylim([0 3])
end

set(gcf,'Position',[684         501        1024         366])