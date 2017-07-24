function P = setupParPriors_MVAR(m)
P.n = 0;
P.params(:,:,1) = randbetween(-2,2,m.n,m.n);
P.params(1,2:3,1) = -32;
P.params(2,1,1) = -32;
P.params(3,2,1) = -32;
P.params(:,:,2) = randbetween(0,0,m.n,m.n).*eye(m.n) + (~eye(m.n).*-32);

P.noisecov      = randbetween(0,0,m.n,m.n).*eye(m.n) + (~eye(m.n).*-32);
