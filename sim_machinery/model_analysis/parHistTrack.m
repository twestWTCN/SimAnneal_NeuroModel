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
parSel = [1 3 5:10];
% Form descriptives
pIndMap = spm_vec(parMu); % in flat form
pIndMap = pIndMap(parSel);
pSigMap = spm_vec(parSigMap); % in flat form
pSigMap = pSigMap(parSel);
parNames = parNames(parSel);
Inds(1,1) = 0; Inds(2,1) = 0;
for multiStart = 1:20
    R.out.dag = sprintf('NPD_STN_GPe_MultiStart_M%.0f',multiStart); % 'All Cross'
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Mfit = varo;
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parHist_' R.out.tag '_' R.out.dag '.mat'])
    parHist = varo;
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\bankSave_' R.out.tag '_' R.out.dag '.mat'])
    bankSave = varo;
    
    parTT = []; r2 = [];
    for i = 1:size(parHist,2)
        parTT(:,i) = spm_vec(parHist(i));
        r2(i) = mean(bankSave{i});
    end
    parT = parTT(pIndMap,:);
    parSig{multiStart} = parTT(pSigMap,:);
    %     parSig = parT(pIndMap,:);
    %     A = 1./parSig;
    %     A = A./sum(A,1);
    A = 1-( parSig{multiStart}.^(1/2));
    %     A = 1/(parSig{multiStart}./parSig{multiStart}(:,1));
    %     A(A<0) = 0
    %     A = 1;
    wParT = A.*parT;
    parWeighted{multiStart} = wParT;
    parMS{multiStart} = parT;
    parConv(:,multiStart) = wParT(:,end);
    Inds(:,multiStart+1) = [Inds(2,multiStart)+1, Inds(2,multiStart) + size(parT,2)];
    R2ms(multiStart) = r2(end);
    R2track{multiStart} = r2;
    Its(multiStart) = size(parT,2);
end
Inds(:,1) = [];
T = [parWeighted{:}];

D = pdist(T','euclidean');
[Y,eigvals] = cmdscale(squareform(D));
format short g
% Check Real Positive
A = [eigvals eigvals/max(abs(eigvals))];

for multiStart = 1:20
    CMDscaled{multiStart} = Y(Inds(1,multiStart):Inds(2,multiStart),:);
end
convMods = find(R2ms>0.01);
convModSel = convMods;
convModSel(convMods>11) = [];
cmap = linspecer(10);
cmap(11:20,:) = repmat([1 0 0],10,1);
i = 0;
figure
subplot(2,2,4)
for multiStart = convMods
    i = i +1;
    %             sc(i) = scatter3(Y(inds,1),Y(inds,2),Y(inds,3),(20.^szvec)*150,cmap(multiStart,:),'.');
    %             hold on
    %             plot3(Y(inds(2:end),1),Y(inds(2:end),2),Y(inds(2:end),3),'color',cmap(multiStart,:))
    
%     sc(i) = scatter(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,3),(15.^R2track{multiStart}(1:3:end))*250,cmap(multiStart,:),'.');
    hold on
    %     p =     plot(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,3),'color',cmap(multiStart,:),'LineWidth',1.25);
    plotVarWidth(CMDscaled{multiStart}(1:3:end,1),CMDscaled{multiStart}(1:3:end,3),2.5.*(15.^R2track{multiStart}(1:3:end)'),cmap(multiStart,:),2)
end
xlabel('Scaling Dimension 1'); ylabel('Scaling Dimension 2');
grid on
title('Multi-start Posteriors on Projected Coordinates')

% legend(sc)

subplot(2,2,3)
% figure
% b = bar(1:numel(pIndMap),mean(parConv,2));
% b.FaceAlpha = 0;
c = bar(1:numel(pIndMap),parConv);
for i = 1:size(c,2)
    c(i).FaceColor = cmap(i,:);
    c(i).EdgeAlpha = 0;
    legnames{i} = sprintf('Model %.0f',i);
end
hold on
b = plot(1:numel(pIndMap),mean(parConv,2),' ko');
b.MarkerFaceColor = 'k';
e = errorbar(1:numel(pIndMap),mean(parConv(:,1:10),2),std(parConv(:,1:10),[],2),'.');
e.Color = [0 0 0];
e.LineWidth = 1.2;
g = gca;
g.XTickLabel = parNames;
g.XTickLabelRotation = -45;
grid on
ylabel('Precision Weighted Posteior Means')
title('Multi-start Posterior Means')

leg = legend(c(1:11),legnames(1:11),'Orientation','vertical');
leg.Box = 'off';
set(leg,'Position',[0.9131 0.1209 0.0803 0.2472]);
set(gcf,'Position',[259         207        1361         615])

for multiStart = convModSel
    % Precision gain
    subplot(2,2,2)
    plot(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}'),2),'color',cmap(multiStart,:),'LineWidth',2);
    hold on
    szvec = R2track{multiStart}; %
    scatter(1:size(parSig{multiStart},2),mean(log10(parSig{multiStart}'),2),(15.^szvec)*100,cmap(multiStart,:),'.');
    % Accuracy gain
    subplot(2,2,1)
    plot(1:size(R2track{multiStart},2),mean(R2track{multiStart}',2),'color',cmap(multiStart,:),'LineWidth',2);
    hold on
    szvec = mean(parSig{multiStart}',2); %
    scatter(1:size(R2track{multiStart},2),mean(R2track{multiStart}',2),(15.^szvec)*100,cmap(multiStart,:),'.');
    
end

subplot(2,2,1)
ylabel('R2');    grid on
xlabel('Iteration');
title('Multi-start Model Convergence')
set(gcf,'Position',[72         296        1831         526])

subplot(2,2,2)
ylabel('Log Precision');    grid on
xlabel('Iteration');
title('Multi-start Parameter Inference')
set(gcf,'Position',[72         296        1831         526])