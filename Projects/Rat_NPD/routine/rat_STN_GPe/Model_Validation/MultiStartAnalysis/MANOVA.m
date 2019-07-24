G = [repmat(1,1,10) repmat(2,1,10)];
X = parConv';
X = parConv(:,1:10)';
X = parConv(:,11:20)';
c = cvpartition(10,'KFold',10);
for i = 1:10
     T2 = T2Hot1(X(~test(c,i),:),0.05,X(test(c,i),:));
     p(i) = T2.p;
end
    
idx = test(c,1);



% T2Hot1(X(1:10,:),0.05)

% for i = 1:8
%     [h p(i)] = ttest(X(:,i));
% end
DepTest1(X(1:10,:),'test','spearman') 


gplotmatrix(X,[],G,[],'+xo');
[d,p,stats] = manova1(X,G);

c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure()
gscatter(c2,c1,G,[],'oxs')
figure()
manovacluster(stats)

mu = mean(X(1:10,:));
sig = cov(X(1:10,:));
R = [mvnrnd(mu,sig,25);  mvnrnd(mu+25,sig+1,25)];
T2Hot1(R,0.05)
