function TimeLockedBetaAnalysis(R,BB,xsimMod,AS)
addpath(R.BBA_path)

PLVcmp = 0;
close all
F(1) = figure; F(2) = figure;
AlignedBursts_analysis(R,1,BB,xsimMod{2}{1},AS{1},PLVcmp,F);
figure(F(1)); title('Full Model');figure(F(2)); title('Full Model')

AlignedBursts_analysis(R,2,BB,xsimMod{1}{1},AS{2},PLVcmp,F);
figure(F(1)); title('No Hyperdirect'); figure(F(2)); title('No Hyperdirect')

AlignedBursts_analysis(R,3,BB,xsimMod{3}{1},AS{3},PLVcmp,F);
figure(F(1)); title('No Cortico-Striatal');figure(F(2)); title('No Cortico-Striatal')

set(F(1),'Position',[374.0000  258.5000  894.5000  620.0000])
set(F(2),'Position',[374.0000  88.0000   894.5000  224.0000])