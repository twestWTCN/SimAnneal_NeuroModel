function BB = compute_BurstThreshold(R,BB,condsel,plotop)
surflag = 0;
% Find the amplitude threshold from the prctile of concatanated data
if surflag == 0
    concatA = [BB.A{condsel}];
    BB.epsAmp = prctile(concatA(R.BB.pairInd(2),:),85,2);
    concatA = [BB.PLV{condsel}];
    BB.epsPLV = prctile(concatA(1,:),90,2);
else
    BORG = load([R.datapathr R.subname{sub} '\ftdata\BetaBursts\BetaBurstAnalysis_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg} '_org'],'BB');
    BB.epsAmp = BORG.BB.epsAmp;
    BB.epsPLV = BORG.BB.epsPLV ;
    clear BORG
end
% Plot Bursting
if plotop == 1
    figure(1)
    plotExampleBurstPLV(R,BB)
    set(gcf,'Position',[680  358  1048  622])
end

save([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'],'BB','R')