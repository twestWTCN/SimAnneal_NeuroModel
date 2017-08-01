function plotDistChange_KS(Rho,nu,xf,psave,R)
% Univariate plots
subplot(1,2,1)
for i = 1:2
    if i == 1
        ls = '--';
        nind = 1;
        
        %{'.params','.noisecov'}
%         M = psave(nind).params; M(M==0) = []; 
        M = psave(nind).A{1}; M(M==0) = []; 
        Ma = M(M>-30);
       
        cmap = linspecer(5);
        X = -5:.1:5;
        for Q = 1:5 %length(Ma)
            p = normpdf(X,Ma(Q),R.SimAn.jitter*R.SimAn.Tm);
%             [p,type,coefs] = pearspdf(X,Ma(Q),R.SimAn.jitter*R.SimAn.Tm,1,3);
            plot(X,p,ls,'color',cmap(Q,:))
            hold on
        end
    else % Copula
        ls = '-';
        r = copularnd('t',Rho,nu,500);
        for Q = 1:5 %size(xf,1)
            x1 = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
            [x1,f] = ksdensity(x1,R.SimAn.pOptRange);
            plot(f,x1,ls,'color',cmap(Q,:))
            hold on
        end
    end
end
    xlabel('\mu')
    ylabel('p(\mu)')
    ylim([0 1]);
    xlim(R.SimAn.pOptBound)
    
    subplot(1,2,2)
    imagesc(Rho)
    set(gca,'YDir','normal')
    title('P Rho')
    set(gca,'XTick',1:size(Rho,1))
    set(gca,'YTick',1:size(Rho,1))
    set(gcf,'Position',[2.5 617 884 383])
    % Multivariate plots
    %
    %         u1 = r(:,1);
    %         v1 = r(:,5);
    %
    %
    %         x1 = ksdensity(xf(1,:),u1,'function','icdf');
    %         y1 = ksdensity(xf(5,:),v1,'function','icdf');
    %
    %         figure;
    %         scatterhist(x1,y1)
    %         set(get(gca,'children'),'marker','.')