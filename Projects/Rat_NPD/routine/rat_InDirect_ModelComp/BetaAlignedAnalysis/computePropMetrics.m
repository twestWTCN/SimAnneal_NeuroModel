function [corrAmp,xcorrAmp,corrPhi,xcorrLAmp,TE,...
    Pv,anTE,peakTau,ZTE] = computePropMetrics(X,Y)

for i = 1:size(X,1)
    for j = 1:size(Y,1)
        XH = hilbert(squeeze(X(i,j,:))'); YH = hilbert(Y(j,:));
        %                         XH_save(:,seg) = abs(XH);
        %                         YH_save(:,seg) = abs(YH);
        ctmp = corrcoef(abs(XH),abs(YH));
        corrAmp(j,i) = ctmp(1,2);
        RP = angle(XH)-angle(YH);
        corrPhi(j,i) =  abs(mean(exp(-1i*RP'),1));
        [c lag] = xcorr(abs(XH),abs(YH),'coeff');
        xcorrAmp(j,i) = max(c);
        xcorrLAmp(j,i) = lag(c==max(c));
        mid = median([XH YH]);
        if nargout>5
        [TE(:,j,i),Pv(:,j,i),anTE(j,i),peakTau(:,j,i),ZTE(:,j,i)] = computeTransferEntropy(squeeze(X(i,j,:)),Y(j,:)',1:30,10);
        end
    end
end
 a = 1;