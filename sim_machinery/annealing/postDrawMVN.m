function par = postDrawMVN(Mfit,pOrg,pIndMap,rep)
disp('Drawing from multivariate normal...')
x1 = mvnrnd(Mfit.Mu',Mfit.Sigma,rep)';

% setup pars from base
clear base
base = repmat(spm_vec(pOrg),1,rep);
for i = 1:rep
    base(pIndMap,i) = x1(:,i);
    par{i} = spm_unvec(base(:,i),pOrg);
end