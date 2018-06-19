function [xint,T] = fx_simulateMVAR(R,x,u,p,m)
% Log Normal Priors
a = m.fn;
P = m;
for i = 1:numel(a)
    y = eval(['m.' a{i}]); ys = eval(['p.' a{i}]);
    Y = y.*exp(ys);
    eval(['P.' a{i} '= Y;'])
end

cfg             = [];
cfg.feedback = 'no';
cfg.triallength = R.IntP.tend;
cfg.ntrials     = 1;
cfg.fsample     = 1/R.IntP.dt;
cfg.nsignal     = m.n;
cfg.method      = 'ar';

a = P.params(:,:,1);
a(a>1) = 1;
cfg.params(:,:,1) = a;

a = P.params(:,:,2);
a(a<-1) = -1;
cfg.params(:,:,2) = a;     

cfg.noisecov      = P.noisecov;

s              = ft_connectivitysimulation(cfg);
xint = s.trial{1};
T = s.time{1};