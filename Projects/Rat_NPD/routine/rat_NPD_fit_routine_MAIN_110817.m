%% RAT DATA- SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% Redo resampling of parameters such that they use the Pind structure to
% map the full P to the opt P.
clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
R = simannealsetup_CSD_app;
rng(8453)
d = R.d; % clock
%% Prepare the data
% % prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\average_rat_lesion.mat'])
% % % Set data as working version

ftdata = ftdata_lesion;
data = ftdata.trial{1}; fsamp = ftdata.fsample;
%Normalise Data
for i = 1:4
    xtmp = data(i,:);
    xtmp = (xtmp-mean(xtmp));%/std(xtmp);
    data(i,:) = xtmp;
end
%compute CSD data features
R.data.fs = ftdata_lesion.fsample;

R.obs.DataOrd = floor(log2(R.data.fs/(2*R.obs.csd.df))); % order of NPD for simulated data
[F_data,meannpd_data] = constructCSDMat(data,R.chloc_name,ftdata.label',fsamp,R.obs.DataOrd,R);

% % [F_data,meannpd_data] = constructNPDMat(data,R.chloc_name,ftdata.label',fsamp,R.obs.DataOrd,R);
save([R.filepathn '\datafeat_npd'],'meannpd_data','F_data')
% %%
load([R.rootn 'data\storage\datafeat_npd']);
R.data.feat_emp = meannpd_data;
R.data.feat_xscale = F_data;
% Plot CSD
csdplotter_220517({meannpd_data},[],F_data,R)
% npdplotter_110717({meannpd_data},[],F_data,R)

%% SIMULATION OF TIME SERIES
% 1 MMC
% 2 STR
% 3 GPE
% 4 STN
% 5 GPI
% 6 THAL

% load priors!
% load([R.rootn 'priors\onerat_bgm_DCM_multipass_split.mat'])
% load([R.rootn 'priors\alt_prior.mat'])
load('C:\Users\twest\Documents\Work\PhD\LitvakProject\Bernadette_DCM\outputs\splitup_BGM\onerat_bgm_DCM_multipass_280617','DCM')
clear u


% Initialise with random priors
% [x,p,m] = setup_new_priors(P,M); % Assign Priors to correct structure

DCM(2) = DCM(1);
%for intrinsic connections
DCM(2).Ep.int{1}.G=DCM(1).Ep.int{1}.G+DCM(1).Ep.int{1}.B;
DCM(2).Ep.int{2}.G=DCM(1).Ep.int{2}.G+DCM(1).Ep.int{2}.B;

%for extrinsic connections
DCM(2).Ep.A{1}=DCM(1).Ep.A{1}+DCM(1).Ep.B{1};
DCM(2).Ep.A{2}=DCM(1).Ep.A{2}+DCM(1). Ep.B{1}; % if only 1 modulatory forward connection (that’s what I use = same modulatory effect thalamus on middle and superficial cortical layers)
% DCM.Ep.A{2}=DCM.Ep.A{2}+DCM.Ep.Bf{2}; % if model contains 2 modulatory forward connections (if there is a second cell in DCM.Ep.Bf)
DCM(2).Ep.A{3}=DCM(1).Ep.A{3}+DCM(1).Ep.B{2};
DCM(2).Ep.A{4}=DCM(1).Ep.A{4}+DCM(1).Ep.B{2}; % if model contains 2 modulatory backward connections (that’s what I use = separate hyperdirect and (in)direct pathway modulations)

% MODEL OFF CONDITION!
DCM = DCM(2);
DCM.Ep.A{1}(6,5) = 0; % Excitatory GPi to Thalamus
%DCM.Ep.A{2}(6,5) = 0;
DCM.Ep.A{1}(1,6) = 0; % Excitatory Thal to M1
% DCM.Ep.A{3}(1,6) = 0; % Inhibitory Thal to M1 
DCM.Ep.A{1}(DCM.Ep.A{1}==0) = rand(size(DCM.Ep.A{1}(DCM.Ep.A{1}==0)));
DCM.Ep.A{3}(DCM.Ep.A{3}==0) = rand(size(DCM.Ep.A{3}(DCM.Ep.A{3}==0)));
% DCM.Ep.A{3}(1,6) = 0; % Inhibitory Thal to M1 

% p.A{3}(6,5) = -32; p.A{4}(6,5) = -32;

% Keep only some fields
p = DCM.Ep;


% Connection strengths
p.C = [-0.0396   -0.1033    0.0025    0.1342   -0.3226    0.0114]'; %zeros(size(p.C,1),1);
p.C_s = repmat(0.5,size(p.C));

% Leadfield
p.obs.LF = zeros(size(R.obs.LF));
p.obs.LF_s = repmat(0.8,size(p.obs.LF));

p.obs.mixing = [-0.2916 -0.0455]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0.4,size(p.obs.mixing));

m = DCM.M;
x = m.x;
A =  p.A; p = rmfield(p,'A');
p.A{1} = A{1}; 
p.A_s{1} = repmat(2.5,size(A{1})); 
p.A{2} = A{3};
p.A_s{2} = repmat(2.5,size(A{3}));

p.S_s = 0.5;
% setup exogenous noise
% m.uset.p = DCM.Ep;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 1e-2; %.*R.InstP.dt;
u = innovate_timeseries(R,m);
u = u./R.IntP.dt;

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.5,size(p.D));

% time constants and gains 
for i = 1:m.m
    if i<5
        p.int{i}.T = zeros(size(p.int{i}.T));
        p.int{i}.T_s = repmat(1,size(p.int{i}.T));
        p.int{i}.G = zeros(size(p.int{i}.G));
        p.int{i}.G_s = repmat(1,size(p.int{i}.G));
        p.int{i}.S = zeros(size(p.int{i}.S));
        p.int{i}.S_s = repmat(1,size(p.int{i}.S));
    else
        p.int{i}.T_s = repmat(0.2,size(p.int{i}.T));
        p.int{i}.G_s = repmat(0.2,size(p.int{i}.G));
        p.int{i}.S_s = repmat(0.2,size(p.int{i}.S));
    end
end
% Correction for STN
% p.int{4}.G = 0.5;
m.n =  size([m.x{:}],2);
m.fxord = DCM.Sname;
m.Bmod = DCM.B;
% set params to zero
% p0 = spm_vec(p);
% p0(p0>-30) = 0;
% 
% p0 = spm_unvec(p0,p);
% p0.B = p.B;

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
tic
%  p = xobs.out.P;
%%%%%%%%%
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\priors\sim_ABC_output_170817_a.mat')
% p = a;
R.objfx.specspec = 'cross'; %%'auto'; % which part of spectra to fit
R.SimAn.jitter = 0.8;
R.out.tag = 'NPD_ABC_autoB4cross';

for i = 1:4
    R.out.tag = [R.out.tag num2str(i)];
%     if i>1
%         p = xobs1.out.P;
%     end
    [xobs1] = SimAn_ABC_110817(x,u,p,m,R);
end
folname = ['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\parfits\' sprintf('%d',[d(1:3)])];
mkdir(folname)
save([folname '\xobs1'],'xobs1');
gif_maker_siman(R)
