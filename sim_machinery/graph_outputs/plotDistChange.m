function [] = plotDistChange(psave,pPrecSave,pSkewSave,stdev,ii)

for i = 1:2
    if i == 1
        ls = '-';
        nind = 1;
    else
        ls = '--';
        nind = ii;
    end

    M = psave(ii-(i-1)).A{1}; M(M==0) = [];
    P = full(pPrecSave{nind}.A{1}); P(P==0) = [];
    S = full(pSkewSave{nind}.A{1}); %S(S==0) = [];
    
    Ma = M(M>-30); 
    P = P(M>-30).*stdev;
    S = S(M>-30);
    
    cmap = linspecer(length(Ma));
    X = -5:.1:5;
    for Q = 1:length(Ma)
        [p,type,coefs] = pearspdf(X,Ma(Q),P(Q)*stdev,S(Q),3);
        plot(X,p,ls,'color',cmap(Q,:))
        hold on
    end
    
end
xlabel('Connection Strength')
ylabel('P(X)')

