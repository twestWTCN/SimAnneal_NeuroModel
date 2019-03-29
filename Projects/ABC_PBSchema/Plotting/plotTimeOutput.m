function plotTimeOutput(xsims_gl,p)
t = xsims_gl.t;
TMS_ind = xsims_gl.TMS_ind;
% tx_iws = xsims_gl.tx_iws;
tx_csn = xsims_gl.tx_csn;
tx_amn = xsims_gl.tx_amn;
tx_emg = xsims_gl.tx_emg;
csn_thresh = xsims_gl.csn_thresh;
amn_thresh = xsims_gl.amn_thresh;

cmap = brewermap(size(tx_amn,2)+2,'YlOrRd');


for i = 1:2
    subplot(2,1,i)
    % TMS_ind(TMS_ind==0) = NaN;
    %   plot(t,tx)
    for cn = 1:size(tx_amn,2)
        pr(cn) = plot(t,tx_csn(:,cn)','color',cmap(cn+1,:));
        hold on
    end
    a(1) = pr(4);
    a(2) = plot(t,repmat(csn_thresh,1,size(t,2)),'g--');
    prc = plot(t,tx_amn','color',[0 0 0]);
    a(3) = prc(1);
    prc = plot(t,repmat(amn_thresh,1,size(t,2)),'b--');
    a(4) = prc(1);
    a(5) = plot(t,tx_emg);
    xlabel('Time (s)')
    ylabel('Amplitude')
    legend(a,'CSN','CSN thresh','AMN','AMN thresh','EMG')
    grid on
    if i == 1
        xlim([0 10])
    else
        xlim([8.5 11])
    end
end