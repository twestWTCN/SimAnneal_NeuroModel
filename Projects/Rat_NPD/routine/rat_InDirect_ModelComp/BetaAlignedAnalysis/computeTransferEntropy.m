function [TE,P,anTE,peakTau,ZTE] = computeTransferEntropy(X,Y,taulist,MCSamp)
% Computes the transfer entropy between X and Y (both directions) by
% scanning over a range of lags (taulist). Statistics are given by MC
% estimation using MCSamps. X must be in matrix dim ch x time (x trial)
% TE(:,1) = Y -> X
% TE(:,2) = X -> Y
if nargin<3
    taulist = 1:10;
end
if nargin<4
    MCSamp = 1000;
end
for tau = taulist
%     [transferEntropy,p,ZTE] = quickTE_TW(X, Y, 'MCnSamples',MCSamp,'delayTrans', tau); 
    meanTE(tau,1) = nan; %nanmean(transferEntropy);
    meanZTE(tau,1) = nan; %nanmean(ZTE);
    meanP(tau,1) =nan; % nanmean(p);
%     [transferEntropy,p,ZTE] = quickTE_TW(Y, X, 'MCnSamples',MCSamp,'delayTrans', tau);
    meanTE(tau,2) = nan; %nanmean(transferEntropy);
    meanZTE(tau,2) = nan; %nanmean(ZTE);
    meanP(tau,2) = nan; %nanmean(p);
end

meanTE(meanTE(:,1)==0,:) = [];
meanZTE(meanZTE(:,1)==0,:) = [];
meanP(meanP(:,1)==0,:) = [];
[a b] = min(meanP);
peakTau = taulist(b);
TE = [meanTE(b(1),1) meanTE(b(2),2)];
ZTE = [meanZTE(b(1),1) meanZTE(b(2),2)];
P = [meanP(b(1),1) meanP(b(2),2)];

anTE = diff(TE)./max(TE);
% anZTE = diff(ZTE)./max(ZTE);