function par = postDrawCopula(R,Mfit,pOrg,pIndMap,pSigMap,rep)
disp('Drawing from copula...')
r = copularnd('t',Mfit.Rho,Mfit.nu,rep);
clear x1
xf = Mfit.xf;
for Q = 1:size(xf,1)
    x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
end
% setup pars from base
clear base
base = repmat(spm_vec(pOrg),1,rep);
for i = 1:rep
    base(pIndMap,i) = x1(:,i);
    base(pSigMap,i) = diag(Mfit.Sigma);
    par{i} = spm_unvec(base(:,i),pOrg);
end

if R.plot.flag == 1
    figure(3)
    clf
    cflag= 1;
    plotProposalDist(R,Mfit,r',pIndMap,cflag)
end