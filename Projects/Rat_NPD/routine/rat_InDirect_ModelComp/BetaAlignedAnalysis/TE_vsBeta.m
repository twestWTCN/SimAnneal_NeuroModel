X = squeeze(TE_mean(1,1,:,8));
Y = squeeze(intpow(2,1,2,:,8));
figure
for p = 1:3
    subplot(2,3,p)
    Y = vertcat(squeeze(intpow(2,1,2,:,:)),nan(1,12));
    Z = vertcat(squeeze(TE_mean(1,p,:,:)),nan(1,12));
    X = repmat([phaseShift],19,1);
    cmap = brewermap(size(Z,2)*2,'Blues');
    
    for i = 1:size(Z,2)
        plot3(X(:,i),Y(:,i),Z(:,i),'LineWidth',2,'Color',cmap(i+12,:));
        hold on
    end
    zlabel('Information Transfer')
    ylabel('Beta Power')
    xlabel('Phase Shift')
    title(coupname{p})
%     if p == 1 || p == 3
%     view([135 35])
%     end
    grid on;% axis equal
    [dum psel] = max(range(Y(1:17)))
   subplot(2,3,p+3)
   plot(Y(:,psel),Z(:,psel),'Color',cmap(psel+12,:),'LineWidth',2)
    ylabel('Information Transfer')
    xlabel('Beta Power')
    grid on; box off;% axis equal
end

set(gcf,'Position',[698         338        1018         597])
