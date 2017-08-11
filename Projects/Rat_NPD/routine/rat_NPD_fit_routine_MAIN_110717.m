%% RAT DATA- SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% Reorganise the parameter shifts - p.C and p.A{2} going HUGE!?!
% ~ Improve the graphic - report the EPS (perct. and while fit), distribution of NMRSEs, 
% ~ also sample copula to show multivariate distribution?
clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
R = simannealsetup_CSD_app;
rng(2042)
d = R.d; % clock
%% Prepare the data
% prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\average_rat_lesion.mat'])
% % % Set data as working version

ftdata = ftdata_lesion;
data = ftdata.trial{1}; fsamp = ftdata.fsample;
%Normalise Data
for i = 1:4
    xtmp = data(i,:);
    xtmp = (xtmp-mean(xtmp))/std(xtmp);
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
load([R.rootn 'priors\onerat_bgm_DCM_multipass_split_finished.mat'])
% load('C:\Users\twest\Documents\Work\PhD\LitvakProject\Bernadette_DCM\outputs\splitup_BGM\onerat_bgm_DCM_multipass_280617','DCM')
clear u
% Initialise with random priors
% [x,p,m] = setup_new_priors(P,M); % Assign Priors to correct structure

DCM(2) = DCM(1);
%for intrinsic connections
DCM(2).Ep.int{1}.G=DCM(1).Ep.int{1}.G+DCM(1).Ep.int{1}.B;
DCM(2).Ep.int{2}.G=DCM(1).Ep.int{2}.G+DCM(1).Ep.int{2}.B;

%for extrinsic connections
DCM(2).Ep.A{1}=DCM(1).Ep.A{1}+DCM(1).Ep.B{1};
DCM(2).Ep.A{2}=DCM(1).Ep.A{2}+DCM(1).Ep.B{1}; % if only 1 modulatory forward connection (that’s what I use = same modulatory effect thalamus on middle and superficial cortical layers)
% DCM.Ep.A{2}=DCM.Ep.A{2}+DCM.Ep.Bf{2}; % if model contains 2 modulatory forward connections (if there is a second cell in DCM.Ep.Bf)
DCM(2).Ep.A{3}=DCM(1).Ep.A{3}+DCM(1).Ep.B{2};
DCM(2).Ep.A{4}=DCM(1).Ep.A{4}+DCM(1).Ep.B{2}; % if model contains 2 modulatory backward connections (that’s what I use = separate hyperdirect and (in)direct pathway modulations)

% MODEL OFF CONDITION!
DCM = DCM(2);
DCM.Ep.A{1}(6,5) = 0; % Excitatory GPe to Thalamus
%DCM.Ep.A{2}(6,5) = 0;
DCM.Ep.A{1}(1,6) = 0; % Excitatory Thal to M1
DCM.Ep.A{3}(1,6) = 0; % Inhibitory Thal to M1 
% DCM.Ep.A{3}(1,6) = 0; % Inhibitory Thal to M1 

% p.A{3}(6,5) = -32; p.A{4}(6,5) = -32;

% Keep only some fields
p = DCM.Ep;
% Connection strengths
p.C = zeros(size(p.C,1),1);
p.C_s = repmat(1,size(p.C));

% Leadfield
p.obs.LF = zeros(1,size(p.C,1));
p.obs.LF_s = repmat(1,size(p.obs.LF));

p.obs.mixing = zeros(1,size(p.C,1));
p.obs.mixing_s = repmat(1,size(p.obs.mixing));

m = DCM.M;
x = m.x;
A =  p.A; p = rmfield(p,'A');
p.A{1} = A{1}; 
p.A_s{1} = repmat(2,size(A{1})); 
p.A{2} = A{3};
p.A_s{2} = repmat(2,size(A{3}));

% setup exogenous noise

% m.uset.p = DCM;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 1; %.*R.IntP.dt;

u = innovate_timeseries(R,m);
u = u./R.IntP.dt;

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.5,size(p.D));

% time constants and gains 
for i = 1:m.m
    p.int{i}.T_s = repmat(1.5,size(p.int{i}.T));
    p.int{i}.G_s = repmat(1.5,size(p.int{i}.G));
    p.int{i}.S_s = repmat(1.5,size(p.int{i}.G));
end

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
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\priors\sim_ABC_output_010817_d.mat')
% p = j;
tic
%  p = xobs.out.P;
for i = 1:4
    R.out.tag = [R.out.tag num2str(i)];
%     if i>1
%         p = xobs1.out.P;
%     end
    [xobs1] = SimAn_ABC_230717(x,u,p,m,R);
end
folname = ['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\parfits\' sprintf('%d',[d(1:3)])];
mkdir(folname)
save([folname '\xobs1'],'xobs1');
gif_maker_siman(R)
% rerun but modulating delays and signal composition
R.SimAn.opPar = {'mix' 'D' 'LF' 'A'};
R.objfx.feattype = 'complex';
R.objfx.specspec = 'cross'; 
R.IntP.intFunc = @stepintegrator_delay;
R.SimAn.saveout = 'xobs2';
R.out.tag = 'crossst';

[xobs1] = SimAn_290617(x,u,p,m,R);
save([folname '\xobs2'],'xobs2');
gif_maker_siman(R)
toc
xsim = stepintegrator(R,x,u,m,xobs2.out.P);


[xsimob,R] = observe_data(xsim,m,xobs2.out.P,R);
[F_sim,meancsd_sim] = constructCSDMat(xsimob,R.chloc_name,R.chloc_name,1/R.IntP.dt,9,R);
figure
csdplotter_220517({meannpd_data},{meancsd_sim},R.frqz,R)
figure
plot(repmat(R.IntP.tvec_obs,size(xsimob,1),1)',xsimob')
legend(R.chsim_name)

figure
stngpe_ob = xsimob(2:3,:);
plot(repmat(R.IntP.tvec_obs,size(stngpe_ob,1),1)',stngpe_ob')
legend(R.chsim_name{2:3})

% % Integrate in time master fx function
 xstore = zeros(6,nt);
% for t = 1:nt5
%     [x xobs] = simAn_master_fx_bgc(x,u,p,m,dt,t);
%     xstore(:,t) = xobs';
%     t
% end
% 
% %% Process Simulated Data
% dobs = [1 4 3 2]; % simulated data matching observed
% % Normalise
% for i = 1:length(dobs)
%     xtmp = xstore(dobs(i),:);
%     xtmp = (xtmp-mean(xtmp))/std(xtmp);
%     xsims(i,:) = xtmp;
% end
% %% Construct CSD
% [F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,ftdata.label',1/dt,10,R.frqz);
% % Plot CSD
% csdplotter([],meancsd_sim,F_sim,R)
