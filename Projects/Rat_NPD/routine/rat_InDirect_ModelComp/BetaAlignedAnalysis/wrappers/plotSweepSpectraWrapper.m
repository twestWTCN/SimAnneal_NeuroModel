function plotSweepSpectraWrapper(R)
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat.mat'],'feat_HD','feat_STR_GPe')
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 6 8 16 18],:);
R.condcmap(6,:) = [0 0 0];
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STR->GPe','150% STR->GPe'};

figure
subplot(1,2,1)
plotSweepSpectra(R.frqz,feat_HD,feat_HD{6},cmap([4 7],:),R.condcmap(1,:),{R.condname{[1 2 3 6]}})
hold on
subplot(1,2,2)
plotSweepSpectra(R.frqz,feat_STR_GPe,feat_STR_GPe{6},cmap([15 18],:),R.condcmap(1,:),{R.condname{[1 4 5 6]}})
set(gcf,'Position',[684         643        1024         224])
