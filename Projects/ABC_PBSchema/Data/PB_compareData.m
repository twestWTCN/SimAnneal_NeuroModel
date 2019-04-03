function r2mean = PB_compareData(R,feat_sim)

for i = 1:size(feat_sim,2)
    yfx = R.data.feat_emp{i}(2,:);
    ffx = feat_sim{i}(2,:);
    %     yfx = yfx - nanmean(yfx);
    %     ffx = ffx - nanmean(ffx);
    r(i) = goodnessOfFit(yfx',ffx','MSE');
%     r2(i) = rsquare(yfx',ffx');
end
r2mean = mean(r([1 3 5 6 7])); %4 5 6 7