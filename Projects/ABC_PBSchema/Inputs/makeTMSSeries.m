function [TMS_onsets,TMS_ind,TMS_win] = makeTMSSeries(t,nPulses,TMS_amp,TMS_dur,TMS_winL,fsamp)
spacing = (fsamp*0.25);
x = spacing + (spacing/4)*randn(1,nPulses);
 x = round(x);
TMS_onsets = cumsum(x); %randi(size(t),1,nPulses); % times of TMS pulses
TMS_onsets = sort(TMS_onsets,'descend');
TMS_onsets(abs(diff(TMS_onsets))<0.1*fsamp) = []; % remove pulses that are too close
TMS_onsets(TMS_onsets>size(t,2)) = [];

TMS_ind = zeros(size(t));

TMS_delay = floor(10.*fsamp/1000);
TMS_delayed_ind = zeros(size(t));
for i = 1:size(TMS_onsets,2)
    TMS_ind(TMS_onsets(i):TMS_onsets(i)+TMS_dur) = repmat(TMS_amp,1,TMS_dur+1);
    TMS_win(TMS_onsets(i):TMS_onsets(i)+TMS_winL) = repmat(1,1,TMS_winL+1);
    TMS_delayed_ind(TMS_onsets(i)+TMS_delay :TMS_onsets(i)+TMS_dur+TMS_delay) = repmat(TMS_amp,1,TMS_dur+1);
end
TMS_delayed_ind = TMS_delayed_ind(1:size(TMS_ind,2));


