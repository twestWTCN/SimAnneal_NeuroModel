function plotTLPhaseSlipProb(TL,cond)
phslip = double(abs(TL.STR_M2_dRP{cond})>0.065);
% phslipN = sum(phslip,1);
% [B,I] = sort(BB.segDur{1});
% phslip = phslip(:,I)
%  subplot(1,2,1)
% imagesc(TL.epochT,1:size(phslip,2),phslip')
% phslipProb =sum(phslip,2)/size(phslip,2);
%  xlim([-300 300])
%  
%  subplot(1,2,2)
plot(TL.epochT(2:end),sum(phslip,2)/size(phslip,2))
xlim([-300 300])
ylabel('Probability of M2/STR Phase Slip')