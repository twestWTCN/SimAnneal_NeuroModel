msrng = [11.385 11.405; 11.82 11.84; 9.91 9.93];
for msz = 1:3
    subplot(3,3,sub2ind([3 3],msz,1))
    imagesc(t,1:p.CSN_n,tx_csn')
    cmap = gray(128); cmap = cmap(end:-1:1,:);
    ylim([1 p.CSN_n])
    colormap(cmap); a(2) = gca;
        grid on
    yyaxis right
    hold on
    plot(t,px_beta); ylim([-pi pi])
    colormap(cmap); a(1) = gca;
    title('CSN Potential')
    
    subplot(3,3,sub2ind([3 3],msz,2))
    % Get Spike Times
    spT = zeros(size(t,2),numel(amn_spT));
    for an = 1:p.AMN_n
        spT(amn_spT{an},an) = 1;
    end
    imagesc(t,1:p.AMN_n,spT')
    colormap(cmap); a(3) = gca;
    ylim([1 p.AMN_n])
    grid on
    title('AMN Units')
    
    subplot(3,3,sub2ind([3 3],msz,3))
    imagesc(t,1,tx_emg')
    colormap(cmap); a(4) = gca;
    grid on
    title('MEP')
    
    linkaxes(a,'x')
    xlim(msrng(msz,:))
end
%hot