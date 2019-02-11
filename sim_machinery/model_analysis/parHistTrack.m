close all
[pInd,parMu,parSig] = parOptInds_110817(R,p,m.m); % in structure form

% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
pSigMap = spm_vec(parSig); % in flat form
Inds = 1;
for multiStart = 1:3
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
    parTT = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
    end
    parT = parTT(pIndMap,:);
    parSig = parTT(pSigMap,:);
%     parSig = parT(pIndMap,:);
    parWeighted = (1./parSig).*parT;
    parMS{multiStart} = parT;
    Inds(multiStart+1) = Inds(multiStart) + (size(parT,2)-1);
end

T = [parMS{:}];

D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));
format short g
% Check Real Positive 
A = [eigvals eigvals/max(abs(eigvals))];
for multiStart = 1:3
    inds = Inds(multiStart:multiStart+1);
    inds = inds(1):inds(2);
    szvec = 1:diff([inds(1) inds(end)])+1;
    scatter3(Y(inds,1),Y(inds,2),Y(inds,3),szvec*10,'.'); 
    hold on
%     plot(Y(inds,1),Y(inds,2)); 
end

