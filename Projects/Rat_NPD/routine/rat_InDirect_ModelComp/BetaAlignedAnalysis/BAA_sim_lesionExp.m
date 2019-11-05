function [R] = BAA_sim_lesionExp(R,simtime)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,10,simtime);
p = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);
 R.Bcond = -1;
cmap = brewermap(18,'Spectral');
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
for i = 1:size(condname,2)
    Pbase_i = Pbase;
    if i~=1
        Pbase_i.A{AIJ{i}(1)}(AIJ{i}(2),AIJ{i}(3)) = -32;
    end
    parfor j = 1:size(condname,2)    % Setup the simulations
        
        Pbase_ij = Pbase_i;
        if j~=1
            Pbase_ij.A{AIJ{j}(1)}(AIJ{j}(2),AIJ{j}(3)) = -32;
        end
        [r2,~,feat_sim] = computeSimData(R,m,uc,Pbase_ij,0);
        powIJ_B(j,i) = sum(feat_sim(1,4,4,3,(R.frqz>14 & R.frqz<30)));
        powIJ_B1(j,i) = sum(feat_sim(1,4,4,3,(R.frqz>14 & R.frqz<21)));
        powIJ_B2(j,i) = sum(feat_sim(1,4,4,3,(R.frqz>21 & R.frqz<30)));
        fitIJ(j,i) = r2;
        disp([i j])
    end
end
 save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','powIJ_B1','powIJ_B2','condname')
% load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','powIJ_B1','powIJ_B2','condname')

% Plot heatmaps
for band = 1
    if band == 1
        X = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100;
        titn = '14-30 Hz';
    elseif band == 2
        X = (powIJ_B1-powIJ_B1(1,1))./(powIJ_B1(1,1)).*100;
        titn = '14-21 Hz';
    elseif band == 3
        X = (powIJ_B2-powIJ_B2(1,1))./(powIJ_B2(1,1)).*100;
        titn = '21-30 Hz';
    end
bcmap = brewermap(128,'RdGy');
subplot(1,2,1)
colormap(bcmap)
igsc = imagesc(X)
set(gca,'YDir','normal');
cb = colorbar;
caxis([-100 100])
title(titn)

% Make Tables

tabdata(:,band) = diag(X);

end

Y = X;
for i = 1:size(X,1)
for j= 1:size(X,2)
    if i~=j
        Y(i,j) = X(i,i)+X(j,j);
    end
end
end

subplot(1,2,2)
colormap(bcmap)
igsc = imagesc(Y-X)
set(gca,'YDir','normal');
caxis([-100 100])

cb = colorbar;

set(gcf,'Position',[1340         124         418         345])


subplot(1,2,2)
colormap(bcmap)
igsc = imagesc(fitIJ)
set(gca,'YDir','normal');
cb = colorbar;

save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','powIJ_B1','powIJ_B2','condname')

% plotSimulatedPRC(PRC,[2 1 3])
% set(gcf,'Position',[711         604        1081         374])