close all
CON = 1;
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'],'BB')

cond = 8
uq = prctile(BB.segDur{cond},85);
lq = prctile(BB.segDur{cond},15);


Xcmb = nchoosek(1:size(BB.epsAmpfull,1), 2);
for i = 1:size(Xcmb,1)
    subplot(3,2,i)
    longRP = BB.segRP{cond}(BB.segDur{cond}>uq,Xcmb(i,1),Xcmb(i,2));
    shortRP = BB.segRP{cond}(BB.segDur{cond}<lq,Xcmb(i,1),Xcmb(i,2));
    
    [pval(i), table] = circ_wwtest(longRP, shortRP);

    % histogram(shortRP,linspace(-pi,pi,12),'Normalization','Probability');
    % hold on
    % histogram(longRP,linspace(-pi,pi,12),'Normalization','Probability');
    %
    polarhistogram(shortRP,linspace(-pi,pi,32),'Normalization','Probability');
    hold on
    polarhistogram(longRP,linspace(-pi,pi,32),'Normalization','Probability');
    
end