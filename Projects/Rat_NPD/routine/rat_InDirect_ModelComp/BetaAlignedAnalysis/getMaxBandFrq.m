function frq_peak = getMaxBandFrq(AS,banddef,Hz)
for i = 1:size(banddef,1)
    for j = 1:size(AS,1)
        [maxpowa fi] = max(AS(j,Hz> banddef(i,1) & Hz< banddef(i,2)));
        frq_peak(i,j) = banddef(1,1)+(fi.*min(diff(Hz)));
    end
end
