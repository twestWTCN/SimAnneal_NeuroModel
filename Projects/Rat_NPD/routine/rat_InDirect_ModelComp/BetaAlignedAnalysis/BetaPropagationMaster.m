function BetaPropagationMaster(R)
addpath(R.BBA_path)
condsel = [5 6 7]; % conditions to be selected
% If doing cond compar
condsel = 1:10;
close all
for CON = 4
    figure
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'],'BB')
     BB.struccmap = linspecer(4);
    TL.struccmap = BB.struccmap;

    ip = 0;
    
    for cond = condsel
        ip = ip + 1;
        TL.periodT = [-250 250];
        %     TL.periodT = [-50 300];
        TL = defineBurstTimeLockEpoch(BB,TL,cond);
        
        figure(1)
        subplot(1,3,ip)
        plotTLTimeEvolutions(TL,cond,'amp',5)
        title(R.condname{cond})
        xlabel('Burst Onset Time (ms)')
        ylim([-20 0])
        xlim([TL.periodT])        
        
        
        figure(2)
        subplot(3,3,ip)
        a = genBoxPlot(TL.onsetT{cond}',BB.struccmap,'horizontal');
        b= gca;
        xlim([TL.periodT])
        b.YTickLabel = R.chsim_name;
        ylabel('Region')
        xlabel('Burst Onset Time (ms)')
        grid on
        title(R.condname{cond})
        
        % Get Onset Stats
        condOnsetStat(:,cond,1) = nanmedian(TL.onsetT{cond}');
        condOnsetStat(:,cond,2) = iqr(TL.onsetT{cond}');
% %         figure(4)
% %         subplot(1,3,ip)
% %         plotTLTimeEvolutions(TL,cond,'BP')
% %         title(R.condname{cond})
% %         xlabel('Burst Onset Time (ms)')
% %         xlim([TL.periodT])
        
                
%         figure(6)
        subplot(3,3,ip+3)
        a = genBoxPlot(TL.maxT{cond}',BB.struccmap,'horizontal');
        b= gca;
        xlim([TL.periodT])
        b.YTickLabel = R.chsim_name;
        ylabel('Region')
        xlabel('Burst Time of 85% Amplitude (ms)')
        grid on
        title(R.condname{cond})
% %         
% %         figure(7)
% %         %     subplot(1,3,ip)
% %         plotTLPhaseSlipProb(TL,cond,'dPhi',ip)
% %         xlim([TL.periodT])
% %         xlabel('Burst Onset Time (ms)')
% %         ylabel('Phase Slip Probability')
% %         xlim([TL.periodT])
% %         grid on
% %         ylim([0 0.5])
% %         
%         figure(8)
        subplot(3,3,ip+6)
        a = genBoxPlot(TL.onsetOffT{cond}',BB.struccmap,'horizontal');
        b= gca;
        xlim([TL.periodT])
        b.YTickLabel = R.chsim_name;
        ylabel('Region')
        xlabel('Time of Theshold Offset(ms)')
        grid on
        title(R.condname{cond})
        
% %         
% %         figure(9)
% %         plotTLPhaseLocking(TL,cond,'dPhi',ip)
% %         xlim([TL.periodT])
% %         xlabel('Burst Onset Time (ms)')
% %         ylabel('Within Burst PLV')
% %         xlim([TL.periodT])
% %         grid on
% %         ylim([0.15 1])
        
    end
    
    set(figure(1),'Position',[680          36        1076         264])
    set(figure(2),'Position',[680          36        1076         942])
%     for i = [3 4 5 6 7 8 9]
%         set(figure(i),'Position',[364         656        1445         322])
%     end
end
% plot(TL.epochT,

% for cond = condsel
%     if isequal(R.condname{cond},'Empirical')
%         T = BB.TEmp;
%         TSw = BB.TSwEmp;
%     else
%         T = BB.T;
%         TSw = BB.TSw;
%     end
%     ip = ip + 1;
%     a(ip) = subplot(1,3,ip)
%     for struc = 1:5
%         if struc<5
%             plot(T(1,1000:end-1000),1*zscore(BB.BPTime{cond}(struc,1000:end-1000)) - struc*10,'color',BB.struccmap(struc,:))
%             %             plot(BB.T(1,1000:end-1000),0.1*(BB.BPTime{cond}(struc,1000:end-1000)) - struc*3,'color',cmap(struc,:))
%         else
%             plot(TSw(1,1000:end-1000),1*zscore(BB.PLV{cond}(1,1000:end-1000)) - struc*10,'color','b')
%             %              plot(BB.TSw(1,1000:end-1000),0.1*(BB.PLV{cond}(1,1000:end-1000)) - struc*3,'color','b')
%         end
%         hold on
%         xlim([120 122])
%     end
% end
% linkaxes(a,'y')
% set(F(5),'Position',[374.0000  88.0000   894.5000  224.0000])
%
%
% F(6) = figure;
% ip = 0;
% for cond = condsel
%     ip = ip + 1;
%     subplot(1,3,ip)
%     for struc = 1:4
%         plot(R.frqz,AS{cond}(struc,:),'color',BB.struccmap(struc,:),'LineWidth',2); hold on
%     end
%     xlabel('Frequency (Hz)'); ylabel('Amplitude')
%     title(R.condname{cond}); grid on
%     ylim([0 35]);
% end
% set(F(6),'Position',[374.0000  88.0000   894.5000  224.0000])
%
% close all
% figure
% a(1) = subplot(5,1,1)
% X(1,:) = wrapToPi(BB.PhiTime{2}(1,:));
% X(2,:) = wrapToPi(BB.PhiTime{2}(2,:));
% X(3,:) = wrapToPi(BB.PhiTime{2}(4,:));
% plot(T,X)
%
% a(2) = subplot(5,1,2)
% XR(1,:) = wrapToPi(BB.PhiTime{2}(1,:)-BB.PhiTime{1}(2,:));
% XR(2,:) = wrapToPi(BB.PhiTime{2}(1,:)-BB.PhiTime{1}(4,:));
% XR(3,:) = wrapToPi(BB.PhiTime{2}(2,:)-BB.PhiTime{1}(4,:));
% plot(T,XR)
% legend({'M1/STR','M1/STN','STR/STN'})
%
% a(3) = subplot(5,1,3)
% XR(1,:) = (BB.PhiTime{2}(1,:)-BB.PhiTime{1}(2,:));
% XR(2,:) = (BB.PhiTime{2}(1,:)-BB.PhiTime{1}(4,:));
% XR(3,:) = (BB.PhiTime{2}(2,:)-BB.PhiTime{1}(4,:));
% XRD = diff(XR,[],2);
% plot(T(2:end),XRD)
%
% a(4) = subplot(5,1,4);
% for struc = 1:5
%     if struc<5
%         plot(T(1,1000:end-1000),1*zscore(BB.BPTime{2}(struc,1000:end-1000)) - struc*10,'color',BB.struccmap(struc,:))
%         %             plot(BB.T(1,1000:end-1000),0.1*(BB.BPTime{cond}(struc,1000:end-1000)) - struc*3,'color',cmap(struc,:))
%     else
%         plot(TSw(1,1000:end-1000),1*zscore(BB.PLV{2}(1,1000:end-1000)) - struc*10,'color','b')
%         %              plot(BB.TSw(1,1000:end-1000),0.1*(BB.PLV{cond}(1,1000:end-1000)) - struc*3,'color','b')
%     end
%     hold on
%     xlim([120 122])
% end
%
% a(5) = subplot(5,1,5)
% plot(T, BB.data{2})
%
% linkaxes(a,'x')