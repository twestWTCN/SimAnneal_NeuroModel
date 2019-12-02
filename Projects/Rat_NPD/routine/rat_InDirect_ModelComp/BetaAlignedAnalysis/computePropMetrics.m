function [corrAmp,xcorrAmp,corrPhi,xcorrLAmp,TE,...
    Pv,anTE,peakTau,ZTE] = computePropMetrics(X,Y)

for i = 1:size(X,1)
    XH = hilbert(X(i,:)); YH = hilbert(Y);
    %                         XH_save(:,seg) = abs(XH);
    %                         YH_save(:,seg) = abs(YH);
    ctmp = corrcoef(abs(XH),abs(YH));
    corrAmp(i) = ctmp(1,2);
    RP = angle(XH)-angle(YH);
    corrPhi(i) =  abs(mean(exp(-1i*RP'),1));
    [c lag] = xcorr(abs(XH),abs(YH),'coeff');
    xcorrAmp(i) = max(c);
    xcorrLAmp(i) = lag(c==max(c));
    mid = median([XH YH]);
    [TE(:,i),Pv(:,i),anTE(i),peakTau(:,i),ZTE(:,i)] = computeTransferEntropy(X(i,:)',Y',1:25,0);
%    TE(:,i) = nan(2,1); Pv = nan(2,1); anTE(i) = nan; peakTau = nan(2,1); ZTE = nan(2,1);
end
