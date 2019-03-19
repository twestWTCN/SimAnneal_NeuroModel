function [out_totP,nspike,nspike_inWin] = spike2psp(inP,thresh,rperiod,ampPSP,epsp,TMS_win)
out_P = nan(size(inP,1),size(inP,2));
for i = 1:size(inP,2)
    in_spT = inP(:,i);
    in_spT = find(in_spT > thresh);
    in_spT(abs(diff(in_spT))<rperiod) = []; % refractory period (50ms)
    nspike = numel(in_spT);
    nspike_inWin = numel(intersect(in_spT,find(TMS_win)));
    P = ampPSP.*convSpikePSP(in_spT,epsp,inP(:,i));
    out_P(:,i) = P;
end
out_totP = sum(out_P,2)';


