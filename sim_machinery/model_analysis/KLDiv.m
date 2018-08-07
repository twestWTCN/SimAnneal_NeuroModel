function [KL DKL] = KLDiv(R,p,m,parBank)
parOptBank = parBank(1:end-1,parBank(end,:)>R.analysis.modEvi.eps);

%% Resample parameters
% Compute indices of optimised parameter
pVec = spm_vec(p);
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form

precInd = parOptIndsPrec_110817(R,p,m.m,2);
precIndMap = spm_vec(precInd);
p = pVec(pIndMap);

R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
xf = zeros(size(pIndMap,1),size(parOptBank,2));
for i = 1:size(pIndMap,1)
    xf(i,:) = parOptBank(pIndMap(i),:); % choose row of parameter values
end

priorVec = spm_vec(R.Mfit.prior);
r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,R.analysis.modEvi.N );
for Q = 1:size(pIndMap,1)
    %     mu_prior(Q) = mean(r(:,Q));
    %     mu_post(Q) = priorVec(pIndMap(Q));
    
    y = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
    [y,f] = ksdensity(y,R.SimAn.pOptRange);
    %     x1
    %     plot(f,x1,'r')
    %     hold on
    x1 = normpdf(f,priorVec(pIndMap(Q)),priorVec(precIndMap(Q)));
    %     plot(f,y,'b')
    
    % Add alpha to get rid of zeros
    y = y + 1e-10;
    x1 = x1 + 1e-10;
    KL(Q) = -sum(y.*log(x1./y));
%     plot(f,x1./y,'g')
end
% KL = sum(KL)
eps1 = priorVec(precIndMap).*eye(size(precIndMap,1));
eps2 = cov(xf');

mu1 = priorVec(pIndMap);
mu2 = mean(xf')';

DKL = 0.5*(trace(eps2\eps1) + ((mu2-mu1)'*eps2)'\(mu2-mu1) - size(mu1,1) + log(det(eps2)/det(eps1)));
