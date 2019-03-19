function [out_spT,nspike,nspike_inWin] = findSpike(inP,thresh,rperiod,TMS_win)
out_spT =[];
for i = 1:size(inP,2)
    in_spT = find(inP(:,i) > thresh);
    in_spT(abs(diff(in_spT))<rperiod) = []; % refractory period (50ms) of incoming burst
    out_spT = union(out_spT,in_spT); % Join together the spike times
%     nspike(i) = numel(in_spT);
    %     nspike_inWin(i) = numel(intersect(in_spT,find(TMS_win)));
end
out_spT(abs(diff(out_spT))<rperiod) = []; % refractory period (50ms) of grouped input (first fires only!)
nspike_inWin = numel(intersect(out_spT,find(TMS_win)));
%     nspike_inWin = size(ismember(out_spT,find(TMS_win)),1);

  nspike= nan;
