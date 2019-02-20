close all; clear
R = simannealsetup_NPD_STN_GPe;
R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',1); % 'All Cross'

modelspec = eval(['@MS_rat_STN_GPe_ModComp_Model' num2str(1)]);
[R,p,m] = modelspec(R);
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;

[pInd,parMu,parSigMap] = parOptInds_110817(R,p,m.m); % in structure form
% [pInd] = createParNameField(R,p,m.m); % in structure form
% Set Par Names
parNames = {'GPe \tau';
    'GPe \gamma';
    'STN \tau';
    'STN \gamma';
    'GPe C';
    'STN C';
    'STN \rightarrow GPe';
    'GPe \rightarrow STN';
    'STN D\rightarrow GPe';
    'GPe D\rightarrow STN';
    };
% Form descriptives
pIndMap = spm_vec(pInd); % in flat form
% pIndMap = pIndMap(5:end);
pSigMap = spm_vec(parSigMap); % in flat form
% pSigMap = pSigMap(5:end);

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
    parSig{multiStart} = parTT(pSigMap,:);
    %     parSig = parT(pIndMap,:);
    %     A = 1./parSig;
    %     A = A./sum(A,1);
    A = 1-( parSig{multiStart}.^(1/2));
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
figure
subplot(1,2,2)
for multiStart = convMods
    i = i +1;
    inds = Inds(multiStart:multiStart+1);
    inds = inds(1):inds(2);
    inds = inds(end-10:end);
    szvec = 1:diff([inds(1) inds(end)])+1;
    %     sc(i) = scatter3(Y(inds,1),Y(inds,2),Y(inds,3),szvec*10,cmap(i,:),'.');
    %     hold on
    %     plot3(Y(inds(2:end),1),Y(inds(2:end),2),Y(inds(2:end),3),'color',cmap(i,:))
    
    sc(i) = scatter(Y(inds,1),Y(inds,2),szvec*10,cmap(i,:),'.');
    hold on
    p =     plot(Y(inds(2:end),1),Y(inds(2:end),2),'color',cmap(i,:));
    legnames{i} = sprintf('Model %.0f',multiStart);
end
xlabel('MDS Dimension 1'); ylabel('MDS Dimension 2');
grid on
title('Multi-start Posteriors on Projected Coordinates')

% legend(sc)

subplot(1,2,1)
% b = bar(1:numel(pIndMap),mean(parConv,2));
% b.FaceAlpha = 0;
c = bar(1:numel(pIndMap),parConv);
for i = 1:size(c,2)
    c(i).FaceColor = cmap(i,:);
end
hold on
b = plot(1:numel(pIndMap),mean(parConv,2),' ko');
b.MarkerFaceColor = 'k';
e = errorbar(1:numel(pIndMap),mean(parConv,2),std(parConv,[],2),'.');
e.Color = [0 0 0];
e.LineWidth = 1.5;
g = gca;
g.XTickLabel = parNames;
g.XTickLabelRotation = -45;
grid on
ylabel('Precision Weighted Posteior Means')
title('Multi-start Posterior Means')

leg = legend(c,legnames,'Orientation','vertical');
leg.Box = 'off';
set(leg,'Position',[0.9131 0.1209 0.0803 0.2472]);
set(gcf,'Position',[259         207        1361         615])

figure
for multiStart = convMods
    plot(mean(log10(parSig{multiStart}'),2),'color',cmap(multiStart,:));
    hold on
end