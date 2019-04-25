function u = KSDensityCVWidth(x,xf,W,range,n,fx)
% Compute using default
[~,~,u] = ksdensity(x,xf,'function',fx,'Weights',W);
% try other bandwidths
uu = logspace(log10(u.*10^range(1)),log10(u.*10^range(2)),n);
v = zeros(size(uu));
% using same partition each time reduces variation 
cp = cvpartition(length(x),'kfold',10);
for j=1:length(uu)
      % compute log likelihood for test data based on training data
      loglik = @(xtr,xte) sum(log(ksdensity(xtr(:,1),xte(:,1),'function',fx,'Weights',xtr(:,3),'width',uu(j))));
      % sum across all train/test partitions
      v(j) = sum(crossval(loglik,[x' xf' W'],'partition',cp));
end 

[~,maxi] = max(v);
u = uu(maxi);
% [f0,x0,u] = ksdensity(x,'width',uu(maxi));
