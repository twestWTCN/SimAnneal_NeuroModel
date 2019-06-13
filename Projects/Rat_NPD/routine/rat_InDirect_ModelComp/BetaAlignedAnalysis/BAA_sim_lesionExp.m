function [R] = BAA_sim_lesionExp(R,MP,simtime)
% Comopute simulations by sweeping across data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
% MP = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);
 R.Bcond = -1;
cmap = brewermap(18,'Spectral');
% PRC.condcmap = cmap([1 4 8 16 4 18],:);
condname = {'Fitted','M2->STR','M2->STN','M2->Thal','STN->GPe','STR-|Gpe','STR->GPi','GPe-|STN','STN->GPi','GPi-|Thal','Thal->M2'};
AIJ = {[],[1 2 1],[1 4 1],[1 6 1],[1 3 4],[2 3 2],[2 5 2],[2 4 3],[1 5 4],[2 6 5],[1 1 6]};

m = MP.m;
P = MP.p;
% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
R.obs.trans.norm = 0;
Pbase = P;
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
        powIJ(j,i) = sum(feat_sim(1,4,4,3,(R.frqz>14 & R.frqz<24)));
        fitIJ(j,i) = r2;
        disp([i j])
    end
end

bcmap = brewermap(128,'RdGy');
subplot(1,2,1)
colormap(bcmap)
igsc = imagesc((powIJ./ powIJ(1,1)).*100)
set(gca,'YDir','normal');
cb = colorbar;

subplot(1,2,2)
colormap(bcmap)
igsc = imagesc(fitIJ)
set(gca,'YDir','normal');
cb = colorbar;

save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ','condname')

% plotSimulatedPRC(PRC,[2 1 3])
% set(gcf,'Position',[711         604        1081         374])