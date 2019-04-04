subplot(3,1,1)
imagesc(t,1:p.CSN_n,tx_csn')
colormap('hot'); a(2) = gca;
yyaxis right
hold on
plot(t,px_beta); ylim([-pi pi])
colormap('hot'); a(1) = gca;
title('CSN Potential')

subplot(3,1,2)
imagesc(t,1:p.AMN_n,tx_amn')
colormap('hot'); a(3) = gca;
title('AMN Potential')

subplot(3,1,3)
imagesc(t,1,tx_emg')
colormap('hot'); a(4) = gca;
title('EMG Potential')

% xlim([101.4 101.45])
linkaxes(a,'x')