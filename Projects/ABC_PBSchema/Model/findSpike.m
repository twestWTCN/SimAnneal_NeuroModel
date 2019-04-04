function [out_spT,nspike,nspike_inWin] = findSpike(inP,thresh,rperiod,TMS_win,t)

% plot(t,inP)
% hold on;

out_spT =[];
for i = 1:size(inP,2)
    in_spT = find(inP(:,i) > thresh);
%     in_spTo = in_spT;
    dL =find(abs(diff(in_spT))<rperiod)+1; % list of differences
    in_spT(dL) = []; % refractory period (50ms) of incoming burst
%     in_spT(abs(diff(in_spT))<rperiod) = []; % refractory period (50ms) of incoming burst
    out_spT{i} = in_spT;
end

out_spT = vertcat(out_spT{:});
nspike_inWin = numel(intersect(out_spT,find(TMS_win)));
nspike= nan;

% tx = zeros(size(inP));
% tx(in_spTo) = 1;
% plot(t,tx.*30)
% 
% tx = zeros(size(inP));
% tx(in_spT) = 1;
% plot(t,tx.*30)
% a =1;
% 
% inP = sum(inP,2);
% in_spT = find(inP > thresh);
% in_spT(abs(diff(in_spT))<rperiod) = []; % refractory period (50ms) of incoming burst
% 
% out_spT = in_spT;
% nspike_inWin = numel(intersect(out_spT,find(TMS_win)));
% nspike= nan;

% for i = 1:size(inP,2)
%     in_spT = find(inP(:,i) > thresh);
%     in_spT(abs(diff(in_spT))<rperiod) = []; % refractory period (50ms) of incoming burst
%     out_spT = union(out_spT,in_spT); % Join together the spike times
% %     nspike(i) = numel(in_spT);
%     %     nspike_inWin(i) = numel(intersect(in_spT,find(TMS_win)));
% end
% out_spT(abs(diff(out_spT))<rperiod) = []; % refractory period (50ms) of grouped input (first fires only!)
% nspike_inWin = numel(intersect(out_spT,find(TMS_win)));
% %     nspike_inWin = size(ismember(out_spT,find(TMS_win)),1);
%
%   nspike= nan;
