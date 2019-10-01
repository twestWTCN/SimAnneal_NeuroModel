function OnsetEvolveAnalysisMaster(R)
addpath(R.BBA_path)
% If doing cond compar
condsel = 1:19;
% close all
a= 0;
for CON = [1 3]
    a = a+1;
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '_bKF.mat'],'BB')
    BB.struccmap = linspecer(4);
    TL.struccmap = BB.struccmap;
    
    condOnsetStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condPeakStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    condTermStat = nan(size(BB.epsAmpfull,1),size(condsel,2),2);
    for cond = condsel
        TL.periodT = [-250 250];
        %     TL.periodT = [-50 300];
        TL = defineBurstTimeLockEpoch(BB,TL,cond);
        
        
        if numel(BB.segInds{cond})>1
            % Get Onset Stats
            condOnsetStat(:,cond,1) = nanmedian(TL.onsetT{cond}');
            condOnsetStat(:,cond,2) = (prctile(TL.onsetT{cond}',66)-prctile(TL.onsetT{cond}',16))./sqrt(size(TL.onsetT{cond},2));
            
            condPeakStat(:,cond,1) = nanmedian(TL.maxT{cond}');
            condPeakStat(:,cond,2) = (prctile(TL.maxT{cond}',66)-prctile(TL.maxT{cond}',16))./sqrt(size(TL.maxT{cond},2));
            
            condTermStat(:,cond,1) = nanmedian(TL.onsetOffT{cond}');
            condTermStat(:,cond,2) = (prctile(TL.onsetOffT{cond}',66)-prctile(TL.onsetOffT{cond}',16))./sqrt(size(TL.onsetOffT{cond},2));
            
            % %         cmblist = nchoosek(1:size(BB.epsAmpfull,1), 2);
            % %         tally = zeros(1,4)
            % %         for p = 1:size(cmblist,1)
            % %             [testR(p),dum,stat] = ranksum(G(:,cmblist(p,1)),G(:,cmblist(p,2)))
            % %             testZ(p) = stat.zval;
            % %         end
            
        end
        
    end
    %     condOnsetStat(isnan(condOnsetStat(:,:,:))) = [];
    subplot(1,2,a)
    betalist = linspace(10,190,19);
    for i = 1:4
        [al aa(i)] = boundedline(betalist,condOnsetStat(i,:,1),[condOnsetStat(i,:,2); condOnsetStat(i,:,2)]')
        al.Color = BB.struccmap(i,:);
        al.LineWidth = 2;
        al.LineStyle = ':';
        aa(i).FaceColor = BB.struccmap(i,:);
        aa(i).FaceAlpha = 0.6;
        
        [bl ba(i)] = boundedline(betalist,condPeakStat(i,:,1),[condPeakStat(i,:,2); condPeakStat(i,:,2)]')
        bl.Color = BB.struccmap(i,:);
        bl.LineWidth = 2;
        bl.LineStyle = '-';
        ba(i).FaceColor = BB.struccmap(i,:);
        ba(i).FaceAlpha = 0.6;
        
        [cl ca(i)] = boundedline(betalist,condTermStat(i,:,1),[condTermStat(i,:,2); condTermStat(i,:,2)]')
        cl.Color = BB.struccmap(i,:);
        cl.LineStyle = '--';
        cl.LineWidth = 2;
        ca(i).FaceColor = BB.struccmap(i,:);
        ca(i).FaceAlpha = 0.6;
        
        %         plot(betalist,condPeakStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2)
        %         hold on
        %         plot(betalist,condOnsetStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2,'LineStyle',':')
        %         plot(betalist,condTermStat(i,:,1),'color',BB.struccmap(i,:),'LineWidth',2,'LineStyle','--')
        
        xlim([10 190]); ylim([-75 130]); grid on
        xlabel('Strength Elliciting %Beta')
        ylabel('Onset Time (ms)')
    end
    set(gca, 'YDir', 'reverse')
    if CON == 1
        set(gca, 'XDir', 'reverse')
    end
end
