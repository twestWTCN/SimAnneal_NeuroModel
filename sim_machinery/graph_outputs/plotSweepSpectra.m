function plotSweepSpectra(Hz,feat,featemp,cmaplims,cmdata,legn)
for c = 1:3
    cmap(:,c) = linspace(cmaplims(1,c),cmaplims(2,c),size(feat,2));
end

for i = 1:size(feat,2)
    if i == 6
        cmapcur = cmdata;
    else
        cmapcur = cmap(i,:);
    end
    a(i) = plot(Hz,squeeze(feat{i}(1,4,4,1,:)),'color',cmapcur,'LineWidth',2);
    hold on
end
a(i+1) = plot(Hz,2e-14.*squeeze(featemp(1,4,4,1,:)),'k:','LineWidth',2);


legend(a([1 6 10]),legn)
xlim([0 48])
xlabel('Frequency (Hz)')
ylabel('Amplitude')
title('Simulated STN Spectra')
box off