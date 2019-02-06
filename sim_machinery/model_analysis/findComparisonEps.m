function R = findComparisonEps(R,r2bank)

prct = 50;
r2bankcat = horzcat(r2bank{:});
R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
% Adjust the acceptance threshold if any models have no rejections
exc = ones(1,numel(R.modcomp.modN));
while any(exc==1)
    r2bankcat = horzcat(r2bank{:});
    R.modcomp.modEvi.epspop = prctile(r2bankcat,prct); % threshold becomes median of model fits
    for modID = 1:numel(R.modcomp.modN)
        exc(modID) = sum(r2bank{modID}>R.modcomp.modEvi.epspop)/size(r2bank{modID},2);
    end
    prct = prct+1;
end