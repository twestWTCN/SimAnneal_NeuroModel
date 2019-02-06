function compStruc = modelEvaluationData(R,modelProbs,modID,compStruc)
% Collate draws
r2rep = [modelProbs.r2rep{:}];
% Remove failed draws
r2rep(isnan(r2rep) | isinf(r2rep)) = [];


compStruc.r2repSave{modID} = (r2rep);
compStruc.KL(modID) = sum(modelProbs.KL(~isnan(modelProbs.KL)));
compStruc.DKL(modID) = sum(modelProbs.DKL);
compStruc.pmod(modID) = sum(r2rep>R.modcomp.modEvi.epspop)/ size(r2rep,2);

list = find([modelProbs.r2rep{:}]>R.modcomp.modEvi.epspop);
if numel(list)>2
    compStruc.parMean{modID} = averageCell(modelProbs.par_rep);
else
    compStruc.parMean{modID} = modelProbs.par_rep{1};
end
