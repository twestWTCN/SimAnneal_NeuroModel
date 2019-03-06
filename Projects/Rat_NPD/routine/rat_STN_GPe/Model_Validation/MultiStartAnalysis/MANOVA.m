G = [repmat(1,1,10) repmat(2,1,10)];
X = parConv';

gplotmatrix(X,[],G,[],'+xo');
[d,p,stats] = manova1(X,G);

c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure()
gscatter(c2,c1,G,[],'oxs')
figure()
manovacluster(stats)