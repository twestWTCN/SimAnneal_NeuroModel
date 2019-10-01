function [pPrec,pSkew] = skewResamp(tbr2,psave,pPrec,pSkew,ii,R)
if ii>5
    n = 5;
else 
    n = 1;
end

for ilag = 1:n
    p1 = spm_vec(psave(ii)); p1(p1==-32) = 0; % p1(p1==0) = NaN;
    p2 = spm_vec(psave(ii-ilag)); p2(p2==-32) = 0; % p2(p2==0) = NaN;
    dTheta(:,ilag) = full(p1-p2);
    dA(:,ilag) = tbr2(ii)-tbr2(ii-ilag);
end

if ii>5
    b = mean(dTheta,2)./mean(dA);
else 
    b = dTheta./dA;
end

skew = -b;
skew = skew.*R.SimAn.dSkew ;
skew = spm_vec(pSkew{ii-1}) + skew; skew(skew>0.7) = 0.4; skew(skew<-0.7) = -0.4;
pSkew{ii} = spm_unvec(skew,psave(ii));

prec = (-abs(b./3))+1;
prec = prec.*R.SimAn.dPrec;
prec = spm_vec(pPrec{ii-1}) + prec; prec(prec<0.25) = 0.4; prec(prec>1.2) = 1.2;
pPrec{ii} = spm_unvec(prec,psave(ii));

% figure
% hist(b)
% 
% for i = -10:2:10
%     r = pearsrnd(0,1,i/10,3,100,1);
%     figure
%     hist(r)
% end

% prec = (-abs(b/3))+1;
% b = linspace(-3,3,25);
% plot(b,prec); shg
end