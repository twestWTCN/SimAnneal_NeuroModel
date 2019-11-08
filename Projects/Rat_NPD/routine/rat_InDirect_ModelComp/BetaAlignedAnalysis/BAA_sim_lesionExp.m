function [R] = BAA_sim_lesionExp(R,simtime,fresh)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,10,simtime);
p = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);
R.Bcond = -1;
cmap = brewermap(128,'RdBu');
ckeypow = linspace(-100,75,128);
ckeyfrq = linspace(-10,10,128);

% PRC.condcmap = cmap([1 4 8 16 4 18],:);
condname = {'Fitted','M2->STR','M2->STN','M2->Thal','STN->GPe','STR-|Gpe','STR->GPi','GPe-|STN','STN->GPi','GPi-|Thal','Thal->M2'};
AIJ = {[],[1 2 1],[1 4 1],[1 6 1],[1 3 4],[2 3 2],[2 5 2],[2 4 3],[1 5 4],[2 6 5],[1 1 6]};

MP.m = m;
MP.p = p;
% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
R.obs.trans.norm = 0;
% R.obs.outstates = find([MP.m.outstates{:}]);
% R.obs.obsstates = 1:4;
Pbase = p;
% powIJ_B = {}; peakIJ_B = {}; freqIJ_B{};
if fresh == 1
    for i = 1:size(condname,2)
        Pbase_i = Pbase;
        if i~=1
            Pbase_i.A{AIJ{i}(1)}(AIJ{i}(2),AIJ{i}(3)) = -32;
        end
        for j = 1:size(condname,2)    % Setup the simulations
            
            Pbase_ij = Pbase_i;
            if j~=1
                Pbase_ij.A{AIJ{j}(1)}(AIJ{j}(2),AIJ{j}(3)) = -32;
            end
            [r2,~,feat_sim] = computeSimData(R,m,uc,Pbase_ij,0);
            [powIJ_B(j,i),peakIJ_B(j,i),freqIJ_B(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,3,:)),[14 30]);
            [powIJ_B1(j,i),peakIJ_B1(j,i),freqIJ_B1(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,3,:)),[14 21]);
            [powIJ_B2(j,i),peakIJ_B2(j,i),freqIJ_B2(j,i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,4,4,3,:)),[21 30]);
            fitIJ(j,i) = r2;
            disp([i j])
        end
    end
    save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
else
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
end
% Plot heatmaps
figure
for band = 1
    titlist = {'14-30 Hz','14-21 Hz','21-30 Hz'}
    if band == 1
        X = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100;
        Y = (freqIJ_B-freqIJ_B(1,1));
    elseif band == 2
        X = (powIJ_B1-powIJ_B1(1,1))./(powIJ_B1(1,1)).*100;
        Y = (freqIJ_B1-freqIJ_B1(1,1));
    elseif band == 3
        X = (powIJ_B2-powIJ_B2(1,1))./(powIJ_B2(1,1)).*100;
        Y = (freqIJ_B2-freqIJ_B2(1,1));
    end
        bcmap = brewermap(128,'RdGy');
        subplot(1,2,1)
        colormap(bcmap)
        igsc = imagesc(X)
        set(gca,'YDir','normal');
        cb = colorbar;
        caxis([-100 100])
%         title(titn)
    
    % Plot as Barplots
    subplot(2,1,1)
    b = bar(diag(X));
    [dum keyind] = min(abs(diag(X)-ckeypow),[],2);
    b.CData = cmap(keyind,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('% Power Change')
    xlim([1.5 11.5])
    subplot(2,1,2)
    b = bar(diag(Y));
    [dum keyind] = min(abs(diag(Y)-ckeyfrq),[],2);
    b.CData = cmap(keyind,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a = gca;
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('Change in Peak Frequency')
    xlim([1.5 11.5])
    % Make Tables
    tabdata(:,band) = diag(X);
    
end
set(gcf,'Position',[448   121   559   592])

Z = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the real (possibly nonlinear effect)
Z(getDiagInd(Z)) = zeros(1,size(Z,1));
X = powIJ_B;
Y = X;
for i = 1:size(X,1)
    for j= 1:size(X,2)
        if i~=j
            Y(i,j) = ((X(i,i)+X(j,j))-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the linear effect
        end
    end
end

figure
colormap(flipud(cmap))
Q = abs(Y)-abs(Z)
Z(isnan(Z)) = 0;
igsc = imagesc(Q)
set(gca,'YDir','normal');
caxis([-125 125]);
a = gca; a.XTick = 1:size(X,1); a.XTickLabel = condname; a.XTickLabelRotation = 45;
a.YTickLabel = condname; a.YTickLabelRotation = 0;
cb = colorbar; axis equal; xlim([1.5 11.5]);ylim([0.5 11.5])

set(gcf,'Position',[1065          84         692         592])


% subplot(1,2,2)
% colormap(cmap)
% igsc = imagesc(fitIJ)
% set(gca,'YDir','normal');
% cb = colorbar;
%
save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','powIJ_B1','powIJ_B2','condname')

% plotSimulatedPRC(PRC,[2 1 3])
% set(gcf,'Position',[711         604        1081         374])