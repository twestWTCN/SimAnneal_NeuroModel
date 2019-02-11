function [Mfit,cflag] = postEstCopula(R,parOptBank,pInd,pIndMap,pOrg)
disp('Forming new copula...')
clear copU xf ilist
% First form kernel density estimates for each optimized
% parameter
clear copU
for i = 1:size(pIndMap,1)
    x = parOptBank(pIndMap(i),:); % choose row of parameter values
    copU(i,:) = ksdensity(x,x,'function','cdf'); % KS density estimate per parameter
    xf(i,:) = x;
end
try
    [Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
    % Save outputs that specify the copula
    Mfit.xf = xf;
    Mfit.ks = copU;
    Mfit.nu = nu;
    Mfit.tbr2 = parOptBank(end,1); % get best fit
    Mfit.Pfit = spm_unvec(mean(parOptBank,2),pOrg);
    Mfit.BPfit = spm_unvec(parOptBank(1:end-1,1),pOrg);
    Mfit.Rho = Rho;
    %     Mfit_hist = Mfit;
    %%% Plot posterior, Rho, and example 2D/3D random draws from copulas
    if R.plot.flag == 1
        figure(3)
        clf
        plotDistChange_KS(Rho,nu,xf,pOrg,pInd,R)
    end
    cflag = 1;
catch
    disp('The estimate of Rho has become rank-deficient.  You may have too few data, or strong dependencies among variables.')
    p = spm_unvec(mean(parOptBank,2),pOrg);
    par = [];
    cflag = 0;
end

Mfit.Mu = mean(parOptBank(pIndMap,:),2);
Mfit.Sigma = cov(parOptBank(pIndMap,:)');


if cflag == 1
    
end