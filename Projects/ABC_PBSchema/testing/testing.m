
[R p m uc] = setup_model_testing(R);
varo = load([R.rootn '\Inputs\betaSig\betaSignal_' num2str(randi(25))]);
m = varo.varo;
varo = load([R.rootn '\Inputs\TMSsig\TMSsignal_' num2str(randi(25))]);
u = varo.varo;

% [xsims dum wflag] = R.IntP.intFx(R,m.x,u,p,m);


[r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,p,0,1);

