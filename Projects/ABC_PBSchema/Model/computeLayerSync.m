function [spSync, spFire, spDelay] = computeLayerSync(windInds,spT)
for i = 1:numel(spT)
    spFire(i) = any(intersect(windInds,spT{i}));
    
    ab = windInds(1) - intersect(windInds,spT{i});
    if any(ab)
        spDelay(i) = ab(1);
    else
        spDelay(i) = NaN;
    end
end

spSync = sum(spFire)/numel(spT);
