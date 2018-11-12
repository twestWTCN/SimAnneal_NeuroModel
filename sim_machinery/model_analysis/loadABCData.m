function [Rout,m,p,parBank] = loadABCData(R)
% Load Options
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\R_' R.out.tag '_' R.out.dag '.mat'])
Rout = varo;
if nargout>1
    % Load Model
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
    m = varo;
end
if nargout>2
    % load modelfit
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
    Rout.Mfit = varo;
    p = Rout.Mfit.BPfit;
end
if nargout>3
    % Load Parbank
    load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
    parBank =  varo;
end