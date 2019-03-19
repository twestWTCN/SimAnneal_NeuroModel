function plotTimeOutput(xsims_gl,p)
t = xsims_gl.t;
TMS_ind = xsims_gl.TMS_ind;
tx_csn = xsims_gl.tx_csn;
tx_amn = xsims_gl.tx_amn;
tx_EMG = xsims_gl.tx_EMG;

cmap = brewermap(size(tx_amn,2)+2,'YlOrRd');


for i = 1:2
    subplot(2,1,i)
    % TMS_ind(TMS_ind==0) = NaN;
    %   plot(t,tx)
    a(1) = plot(t,TMS_ind,'LineWidth',2);
    hold on
    a(2) = plot(t,repmat(p.SP_eps(2),1,size(t,2)),'k--');
    for cn = 1:size(tx_amn,2)
        pr(cn) = plot(t,tx_csn(:,cn)','color',cmap(cn+1,:));
    end
    a(3) = pr(4);
    a(4) = plot(t,tx_amn(:,4),'color',[0 0 0]);
    a(5) = plot(t,tx_EMG);
    xlabel('Time (s)')
    ylabel('Amplitude')
    legend(a,'TMS','CSN Threshold','CSN','AMN','EMG')
    grid on
    if i == 1
        xlim([20 30])
    else
        xlim([23.5 25])
    end
end