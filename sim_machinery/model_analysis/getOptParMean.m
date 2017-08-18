function optP = getOptParMean(m,p,R,d)
%% Resample parameters
load([R.rootn 'outputs\' R.out.tag '\parBank_' R.out.tag '_' d '.mat'])
parOptBank = varo;
% figure
% hist(parOptBank(end,:),[-1:.1:1]); xlim([-1 1])
eps = R.analysis.modEvi.eps;
N = R.analysis.modEvi.N;

% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
xf = zeros(size(pIndMap,1),size(parOptBank,2));
for i = 1:size(pIndMap,1)
    x = parOptBank(pIndMap(i),:); % choose row of parameter values
    xf(i,:) = x;
end

disp('Drawing from copula...')
r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,N);
clear x1
for Q = 1:size(xf,1)
    x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
end
% setup pars from base
clear base
base = repmat(spm_vec(p),1,N);
for i = 1:N
    base(pIndMap,i) = x1(:,i);
end
optP = spm_unvec(mean(base,2),p);