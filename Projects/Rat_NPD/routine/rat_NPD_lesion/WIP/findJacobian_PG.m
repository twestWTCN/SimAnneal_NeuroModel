%function fndJacobian(x,m,p,R)

%% Simulate New Data
% u = innovate_timeseries(R,m);
% u{1} = u{1}.*sqrt(R.IntP.dt);
% xsims = spm_fx_compile_findJacobian(R,m.x,u,p,m);

R.IntP.nt = R.IntP.buffer + R.IntP.buffer;
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);

u = {zeros(size(u{1}))};

fx = @(x)spm_fx_compile_findJacobian(x,R,u,p,m);
options.PlotFcns = @optimplotx
fminsearch(fx,spm_vec(m.x))