function plotSweepSpectraWrapper(R)
close all
load([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat_F1.mat'],'feat_HD','feat_STR_GPe')
cmap = brewermap(30,'Spectral');
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(30,'Reds');
% cmap1(22,:) = [0 0 0];
cmap2 = brewermap(30,'Blues');
% cmap2(28,:) = [0 0 0];

% condrange = 

figure
subplot(1,2,1)
plotSweepSpectra(R.frqz,feat_HD,feat_HD{6},cmap1,{R.condname{[2 1 3]}},[1 22 25],1:3:30)
hold on
subplot(1,2,2)
x = [3:2:30];x(end-2) = 27; x(end-2) = 28; x(end-1) = 29;
plotSweepSpectra(R.frqz,feat_STR_GPe,feat_STR_GPe{6},cmap2,{R.condname{[5 4 6]}},[3 28 29],x)
set(gcf,'Position',[684         501        1024         366])


ck_1 = logspace(-1,0.6,30);
ck_2 = logspace(-1,0.14,30);

% cmap = brewermap(30,'Reds');


figure
for ck = 1:30
bpow(ck,1) = 1e7.*max(feat_HD{ck}(1,4,4,3,R.frqz>14 & R.frqz<30));
bpow(ck,2) = 1e7.*max(feat_STR_GPe{ck}(1,4,4,3,R.frqz>14 & R.frqz<30));
[a b] = max(feat_HD{ck}(1,4,4,3,R.frqz>14 & R.frqz<30));
fpow(ck,1) = 14+ R.frqz(b);

[a b] = max(feat_STR_GPe{ck}(1,4,4,3,R.frqz>14 & R.frqz<30));
fpow(ck,2) = 14+ R.frqz(b);
end

subplot(1,2,1)
br = plot(log10(ck_1(1,:)),log10(bpow(:,1)),'k-');
hold on
Sr = scatter(log10(ck_1(1,:)),log10(bpow(:,1)),50,cmap1,'filled');
ylabel('Peak Amplitude (V)')
grid on

yyaxis right
br = plot(log10(ck_1(1,:)),(fpow(:,1)),':');
hold on
Sr = scatter(log10(ck_1(1,:)),(fpow(:,1)),50,cmap1,'filled');
Sr.Marker = 'diamond';
ylabel('Peak Frequency (Hz)')
xlabel('log_{10} % Connection Strength')
title('Modulating M2 -> STN')

subplot(1,2,2)
br = plot(log10(ck_2(1,:)),log10(bpow(:,2)),'k-');
hold on
Sr = scatter(log10(ck_2(1,:)),log10(bpow(:,2)),50,cmap2,'filled');
ylabel('Peak Amplitude (V)')

yyaxis right
br = plot(log10(ck_2(1,:)),(fpow(:,2)),':');
hold on
Sr = scatter(log10(ck_2(1,:)),(fpow(:,2)),50,cmap2,'filled');
Sr.Marker = 'diamond';
ylabel('Peak Frequency (Hz)')
xlabel('log_{10} % Connection Strength')
title('Modulating STN -| GPe')
grid on

set(gcf,'Position',[684         501        1024         366])