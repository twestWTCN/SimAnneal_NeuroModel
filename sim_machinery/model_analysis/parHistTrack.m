close all
[pInd,parMu,parSig] = parOptInds_110817(R,p,m.m); % in structure form
% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
pIndMap = [3 6 39 41 45 47];
Inds = 1;
for multiStart = 1:6
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
    parT = [];
    for i = 1:size(parHist,2)
        parT(:,i) = spm_vec(parHist(i));
    end
    parT = parT(pIndMap,:);
    parMS{multiStart} = parT;
    Inds(multiStart+1) = Inds(multiStart) + (size(parT,2)-1);
end

T = [parMS{:}];

D = pdist(T','euclidean');
squareform(D)
[Y,eigvals] = cmdscale(D);
format short g
[eigvals eigvals/max(abs(eigvals))]
for multiStart = 1:6
    inds = Inds(multiStart:multiStart+1);
    inds = inds(1):inds(2);
    szvec = 1:diff([inds(1) inds(end)])+1;
    scatter(Y(inds,1),Y(inds,2),szvec*10,'.'); 
    hold on
%     plot(Y(inds,1),Y(inds,2)); 
end

