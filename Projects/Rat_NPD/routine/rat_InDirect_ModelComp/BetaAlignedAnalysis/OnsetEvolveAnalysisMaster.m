function OnsetEvolveAnalysisMaster(R)
addpath(R.BBA_path)
% If doing cond compar
condsel = 1:10;
% close all
for CON = 1:4
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'],'BB')
    BB.struccmap = linspecer(4);
    TL.struccmap = BB.struccmap;
    
    condOnsetStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    
    for cond = condsel
        TL.periodT = [-250 250];
        %     TL.periodT = [-50 300];
        TL = defineBurstTimeLockEpoch(BB,TL,cond);
        
        
        if numel(BB.segInds{cond})>1
            % Get Onset Stats
            condOnsetStat(:,cond,1) = nanmedian(TL.onsetT{cond}');
            condOnsetStat(:,cond,2) = (prctile(TL.onsetT{cond}',66)-prctile(TL.onsetT{cond}',16))./sqrt(size(TL.onsetT{cond},2));
            
            % %         cmblist = nchoosek(1:size(BB.epsAmpfull,1), 2);
            % %         tally = zeros(1,4)
            % %         for p = 1:size(cmblist,1)
            % %             [testR(p),dum,stat] = ranksum(G(:,cmblist(p,1)),G(:,cmblist(p,2)))
            % %             testZ(p) = stat.zval;
            % %         end
            
        end
        
    end
    %     condOnsetStat(isnan(condOnsetStat(:,:,:))) = [];
    subplot(4,1,CON)
    for i = 1:4
        [bl ba(i)] = boundedline(log10(BB.condlist),condOnsetStat(i,:,1),condOnsetStat(i,:,2))
        bl.Color = BB.struccmap(i,:);
        bl.LineWidth = 2;
        ba(i).FaceColor = BB.struccmap(i,:);
        ba(i).FaceAlpha = 0.6;
        hold on
        xlim([-1 1])
        xlabel('log Connection Strength')
        ylabel('Onset Time (ms)')
    end
    if CON == 4
        legend(ba,R.chsim_name)
    end
end
