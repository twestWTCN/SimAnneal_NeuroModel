function PB_PlotOutputs(R,feat_emp,feat_sim,bestn)
if isempty(feat_emp)
    for i = 1:size(feat_sim,2)
        feat_emp{1}{i} = zeros(size(feat_sim{1}{i}));
    end
end
sflag = 1;
if isempty(feat_sim)
    for i = 1:size(feat_emp,2)
        sflag = 0;
        feat_sim{1}{i} = zeros(size(feat_emp{i}));
    end
end

plotTMS_MEPStats(feat_emp{1},[0 0 0],1)
if sflag ==1
    for L = 1:size(feat_sim,2)
        if L == 1
            lwid = 2;
        else
            lwid = 0.5;
        end
        
        plotTMS_MEPStats(feat_sim{L},[1 0 0],lwid)
    end
end
 set(gcf,'Position',[1210          23         546         923])