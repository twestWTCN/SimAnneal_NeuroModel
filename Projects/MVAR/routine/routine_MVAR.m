%% MVAR - SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% CURRENTLY NO AUTOSPECTRA??!
clear ; close all

%%ft_postamble_provenance HAS BEEN ALTERED
%% Create 'data features' to be fitted
R = simannealsetup_MVAR;

M = setupParVals_MVAR; % values (model)
M.m = 1; % 1 "compartment"
P = setupParPriors_MVAR_orig(M); % Prior means)
x = [];

[xint T] = fx_simulateMVAR(x,M,P,R);
[F meannpd] = R.obs.transFx(xint,R.chloc_name,R.chloc_name,1/R.IntP.dt,7,R);

R.data.feat_emp = meannpd;
R.data.feat_xscale = F;
R.plot.outFeatFx({R.data.feat_emp},{},R.data.feat_xscale,R,1)

clear xint T P


%% Try and infer the parameters of the generative model using simulated annealing
x = [];
P = setupParPriors_MVAR(M);
% [xobs1] = SimAn_100717(x,[],P,M,R);
[R] = SimAn_ABC_230717(x,[],P,M,R)
% gif_maker_siman(R)