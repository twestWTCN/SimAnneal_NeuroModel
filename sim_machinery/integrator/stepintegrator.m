function xstore = stepintegrator(R,x,u,m,p)
xstore = zeros(m.n,R.IntP.nt);
for t = 1:R.IntP.nt
    [xint] = spm_fx_gen2(x,u(t,:),p,m,R.IntP.dt);
    xstore(:,t) = xint;
    x = spm_unvec(xint,m.x);
    t;
end