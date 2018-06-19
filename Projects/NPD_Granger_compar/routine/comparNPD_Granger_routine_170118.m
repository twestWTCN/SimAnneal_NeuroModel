%% RAT DATA- SIM ANNEAL PROJECT
%%%%%%%%%%%%%%%%%%
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
rng(623)

%% Set Parameters of the Routine
R = simannealsetup_comparNPD_Granger_170118;
R.d = clock; % clock

%% Prepare Model
m.m = 2; % # of sources
m.x = {[0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0]}; % Initial states
m.Gint = [8 8] ; % Priors
m.Tint = [4 4];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
m.dipfit.model(1).source = 'CMC';
m.dipfit.model(2).source = 'CMC';

% Precompute xinds to make things easier with indexing
% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
        xinds(i,1) = 1;
        xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;

% setup exogenous noise
% m.uset.p = DCM.Ep;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 1e-2; %.*R.InstP.dt;
rng(8453)
u = innovate_timeseries(R,m);
u = u./R.IntP.dt;


%% Prepare Priors
% 1 X1
% 2 X2


% Excitatory connections
A = repmat(-32,m.m,m.m);
% A(2,1) = -12;
A(1,2) = -8;

p.A{1} = A;
A_s = repmat(0.2,size(A));
p.A_s{1} = A_s;

% Inhbitory connections
A = repmat(-32,m.m,m.m);
% A(2,1) = 6;
p.A{2} = A;
A_s = repmat(0.2,size(A));
p.A_s{2} = A_s;

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.2,size(p.C));

% Leadfield
p.obs.LF = zeros(size(R.obs.LF));
p.obs.LF_s = repmat(0.2,size(p.obs.LF));

p.obs.mixing = [0 0]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0.2,size(p.obs.mixing));

% Delays
p.D = repmat(2,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.01,size(p.D));

% Sigmoid transfer for connections
p.S = 0;
p.S_s = 0.01;

% time constants and gains
for i = 1:m.m
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(1,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(1,size(p.int{i}.G));
    p.int{i}.S = zeros(1);
    p.int{i}.S_s = repmat(1,size(p.int{i}.S));
end

%%%%%%%%%
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\priors\sim_ABC_output_170817_a.mat')
% p = a;
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\cross_fit_workspace.mat')
% R = R_out;
R.obs.states = [8 16];
R.objfx.specspec = 'cross'; %%'auto'; % which part of spectra to fit
R.SimAn.jitter = 0.5;
m.uset.p.scale = 0.05;
R.out.tag = 'comparNPD_Granger';

rng(124)
p = resampleParameters_240717(R,p,1,m.m); % Draw from prior
xsims = R.IntP.intFx(R,m.x,u,p,m); %spm_fx_compile_180817
xsims_gl = R.obs.obsFx(xsims,m,p,R);

figure(1)
tvec_obs = R.IntP.tvec;
tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
R.IntP.tvec_obs = tvec_obs;
subplot(2,1,1)
plot(repmat(R.IntP.tvec_obs,size(xsims_gl,1),1)',xsims_gl');
xlabel('Time (s)'); ylabel('Amplitude')
subplot(2,1,2)
plot(repmat(R.IntP.tvec_obs,size(xsims_gl,1),1)',xsims_gl'); xlim([5 6])
xlabel('Time (s)'); ylabel('Amplitude')
legend(R.chsim_name)
drawnow;shg

% normnoise(1,:) = sin(2*pi*25.*tvec_obs)+ (rand(1,length(xsims_gl))/2);
% normnoise(2,:) = sin(2*pi*25.*tvec_obs + pi/2)+ (rand(1,length(xsims_gl))/2);
cfg             = [];
cfg.ntrials     = 10;
cfg.triallength = tvec_obs(end);
cfg.fsample     = 1/R.IntP.dt;
cfg.nsignal     = 2;
cfg.method      = 'ar';
cfg.params(:,:,1) = [
    1   0;
    0   1];

cfg.params(:,:,2) = [
    -0.5  0 ;
    0   -0.8 ] ;
cfg.params(:,:,3) = [
    0  0 ;
    -1  0] ;
cfg.noisecov      = [
    0.001    0   ;
    0    0.001 ] ;
data              = ft_connectivitysimulation(cfg);
tvec = data.time{1};
normnoise = data.trial{1};

[F{1},feat_sim{1}] = constructNPDMat_180118(xsims_gl,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R,normnoise);

%         figure(2)
%         clf
%         fx = R.plot.outFeatFx;
%         fx({},feat_sim,F,R,1,[])
%         drawnow; shg

[F{2},feat_sim{2}] = construct_spectral_ARGrangeMat(xsims_gl,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R,normnoise);
%         figure(3)
%         clf
%         fx = R.plot.outFeatFx;
%         fx({},feat_sim,F,R,1,'AR Spectral Granger')
%         drawnow; shg
% [F{3},feat_sim{3}] = construct_spectral_NPGrangeMat(xsims_gl,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R,normnoise);
%         figure(4)
%         clf
%         fx = R.plot.outFeatFx;
%         fx({},feat_sim,F,R,1,'NP Spectral Granger')
%         drawnow; shg
% [F{3} feat_sim{3}] = construct_spectral_NP_DTFMat(xsims_gl,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R,normnoise);
%         figure(5)
%         clf
%         fx = R.plot.outFeatFx;
%         fx({},feat_sim,F,R,1,'NP Spectral DTF')
%         drawnow; shg
figure(6)
clf
fx = R.plot.outFeatFx;
comparplotter_180117({},feat_sim,F,R,1,'Metric',{'NPD','AR Granger'})
drawnow; shg