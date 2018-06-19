function P = setupParPriors_MVAR(m)
P.n = 0;
P.params(:,:,1) = randbetween(-5,5,m.n,m.n);
P.params(1,2:3,1) = -32;
% P.params(2,1,1) = -32;
P.params(3,2,1) = -32;
P.params(:,:,2) = randbetween(3,3,m.n,m.n).*eye(m.n) + (~eye(m.n).*-32);
P.params_s = repmat(2,size(P.params));
P.noisecov      = randbetween(3,3,m.n,m.n).*eye(m.n) + (~eye(m.n).*-32);
P.noisecov_s = repmat(2,size(P.params));

