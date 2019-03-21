subplot(3,1,1)
title('CSN Potential')
imagesc(t,1:CSN_n,tx_csn')
colormap('hot'); a(2) = gca;
yyaxis right
hold on
plot(t,px_beta); ylim([-pi pi])
colormap('hot'); a(1) = gca;

subplot(3,1,2)
title('AMN Potential')
imagesc(t,1:AMN_n,tx_amn')
colormap('hot'); a(3) = gca;

subplot(3,1,3)
title('EMG Potential')
imagesc(t,1,tx_emg')
colormap('hot'); a(4) = gca;

% xlim([101.4 101.45])
linkaxes(a,'x')