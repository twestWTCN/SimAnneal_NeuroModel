close all
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
Mfit = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;

varo = load([R.rootn '\Inputs\TMSsig\TMSsignal_' num2str(randi(25))]);
u = varo.varo;
varo = load([R.rootn '\Inputs\betaSig\betaSignal_' num2str(randi(25))]);
m = varo.varo;


clear varo

BPfit = Mfit.BPfit;
% [r2,pnew,feat_sim,xsims,xsims_gl] = computeSimData120319(R,m,u,BPfit,0,1);
R = PB_scheme_setup();

BPfit.SCRate(1) = 0; BPfit.SCRate(2) =0; BPfit.SCRate(3) = 0;
 [xsims,dum,wflag] = PB_schema_simulate_REV(R,[],u,BPfit,m);