function plotSweepSpectra(Hz,feat,featemp,cmap,legn,legsel,condsel)

for i = condsel
    a(i) = plot(Hz,1e7.*squeeze(feat{i}(1,4,4,1,:)),'color',cmap(i,:),'LineWidth',2);
    hold on
end
% a(i+1) = plot(Hz,2e-14.*squeeze(featemp(1,4,4,1,:)),'k:','LineWidth',2);

a(legsel(2)).LineStyle = '--';
legend(a(legsel),legn)
xlim([0 48])
xlabel('Frequency (Hz)')
ylabel('Amplitude (uV Hz^-1)')
title('Simulated STN Spectra')
box off
set(gca, 'YScale', 'log'); %, 'XScale', 'log')
grid on