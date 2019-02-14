close all; clear
R = simannealsetup_NPD_STN_GPe;
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'

modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
[R,p,m] = modelspec(R);
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;

[pInd,parMu,parSig] = parOptInds_110817(R,p,m.m); % in structure form

% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
pIndMap = pIndMap(5:end);
pSigMap = spm_vec(parSig); % in flat form
pSigMap = pSigMap(5:end);

Inds = 1;
for multiStart = 1:10
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Mfit = varo;
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
    parTT = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
    end
    parT = parTT(pIndMap,:);
    parSig = parTT(pSigMap,:);
    %     parSig = parT(pIndMap,:);
    A = 1./parSig;
    A = A./sum(A,1);
    wParT = A.*parT;
    parWeighted{multiStart} = wParT;
    parMS{multiStart} = parT;
    parConv(:,multiStart) = wParT(:,end);
    Inds(multiStart+1) = Inds(multiStart) + (size(parT,2)-1);
    R2ms(multiStart) = Mfit.tbr2;
end

T = [parWeighted{:}];

D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));
format short g
% Check Real Positive
A = [eigvals eigvals/max(abs(eigvals))];

convMods = find(R2ms>-1);
cmap = linspecer(10);
i = 0;
for multiStart = convMods
    i = i +1;
    inds = Inds(multiStart:multiStart+1);
    inds = inds(1):inds(2);
    szvec = 1:diff([inds(1) inds(end)])+1;
%     sc(i) = scatter3(Y(inds,1),Y(inds,2),Y(inds,3),szvec*10,cmap(i,:),'.');
%     hold on
%     plot3(Y(inds(2:end),1),Y(inds(2:end),2),Y(inds(2:end),3),'color',cmap(i,:))
    
    sc(i) = scatter(Y(inds,1),Y(inds,2),szvec*10,cmap(i,:),'.');
    hold on
    p =     plot(Y(inds(2:end),1),Y(inds(2:end),2),'color',cmap(i,:));
end
legend(sc)

figure
bar(pIndMap,mean(parConv,2))
hold on
errorbar(pIndMap,mean(parConv,2),std(parConv,[],2),'.')
