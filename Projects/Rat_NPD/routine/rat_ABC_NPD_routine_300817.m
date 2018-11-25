%% RAT DATA- SIM ANNEAL PROJECT
%%%%%%%%%%%%%%%%%%
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
rng(231231)

%% Set Parameters of the Routine
R = simannealsetup_NPD_300817;
R.d = clock; % clock

%% Prepare the data
% prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\NPD_paper.mat']);
X = squeeze(meannpd_data(:,:,:,:,2,:));
meannpd_data = squeeze(mean(X,5));
for i = 1:size(meannpd_data,1)
    for j = 1:size(meannpd_data,2)
        for k = 1:size(meannpd_data,3)
            Pxy = squeeze(meannpd_data(i,j,k,:));
            Pxy = interp1(F_data,Pxy,R.frqz);
            Pxy = Pxy.*tukeywin(length(Pxy),0.25)';
            tmp(i,j,k,:) = Pxy;
        end
    end
end
% Set data as working version
R.data.feat_emp = tmp;
R.data.feat_xscale = R.frqz;

% Plot CSD
if strcmp('CSD',R.data.datatype)
    csdplotter_220517({meannpd_data},[],F_data,R)
elseif strcmp('NPD',R.data.datatype)
    npdplotter_110717({meannpd_data},[],F_data,R)
end

%% Prepare Model
m.m = 6; % # of sources
m.x = {[0 0 0 0 0 0 0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]}; % Initial states
m.Gint = [14 1 1 1 1 1];
m.Tint = [4 1 1 1 1 1];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
m.dipfit.model(1).source = 'MMC';
m.dipfit.model(2).source = 'STR';
m.dipfit.model(3).source = 'GPE';
m.dipfit.model(4).source = 'STN';
m.dipfit.model(5).source = 'GPI';
m.dipfit.model(6).source = 'THAL';

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
u = innovate_timeseries(R,m);
u = u./R.IntP.dt;


%% Prepare Priors
% 1 MMC
% 2 STR
% 3 GPE
% 4 STN
% 5 GPI
% 6 THAL

% Excitatory connections
A = repmat(-32,m.m,m.m);
A(2,1) = 0; % M1 to STR
A(4,1) = 0.5; % M1 to STN
A(3,4) = 0.5; % STN to GPe
A(5,4) = 0; % STN to GPi
A(1,6) = 0; % THAL to M1
p.A{1} = A;
A_s = repmat(1,size(A));
p.A_s{1} = A_s;

% Inhbitory connections
A = repmat(-32,m.m,m.m);
A(3,2) = 0.5; % STR to GPe
A(5,2) = 0; % STR to GPi
A(4,3) = 0.5; % GPe to STN
A(5,3) = 0; % GPe to GPi
A(2,3) = 0; % GPe to STR
A(6,5) = 0; % GPi to THAL
p.A{2} = A;
A_s = repmat(1,size(A));
p.A_s{2} = A_s;

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.5,size(p.C));

% Leadfield
p.obs.LF = zeros(size(R.obs.LF)).*0.8;
p.obs.LF_s = repmat(0.2,size(p.obs.LF));

p.obs.mixing = [-0.2916 -0.0455]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0.2,size(p.obs.mixing));

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.8,size(p.D));

% Sigmoid transfer for connections
p.S = 0;
p.S_s = 1;

% time constants and gains
for i = 1:m.m
    if i<5
        p.int{i}.T = zeros(1,m.Tint(i));
        p.int{i}.T_s = repmat(0.5,size(p.int{i}.T));
        p.int{i}.G = zeros(1,m.Gint(i));
        p.int{i}.G_s = repmat(0.5,size(p.int{i}.G));
        p.int{i}.S = zeros(1);
        p.int{i}.S_s = repmat(0.5,size(p.int{i}.S));
    else
        p.int{i}.T = zeros(1,m.Tint(i));
        p.int{i}.T_s = repmat(0.5,1,m.Tint(i));
        p.int{i}.G = zeros(1,m.Gint(i));
        p.int{i}.G_s = repmat(0.5,1,m.Gint(i));
        p.int{i}.S = zeros(1);
        p.int{i}.S_s = repmat(0.5,1,1);
    end
end

%%%%%%%%%
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\priors\sim_ABC_output_170817_a.mat')
% p = a;
% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\cross_fit_workspace.mat')
% R = R_out;
R.objfx.specspec = 'auto'; %%'auto'; % which part of spectra to fit
R.plot.save = 'True';
R.SimAn.jitter = 0.4;
m.uset.p.scale = 0.05;

for i = 1:4
    %     R.out.tag = [R.out.tag num2str(i)];
    if i>1
        p = R_out.Mfit.Pfit;
    end
    [xobs1] = SimAn_ABC_110817(m.x,u,p,m,R);
end
% folname = ['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\parfits\' sprintf('%d',[d(1:3)])];
% mkdir(folname)
% save([folname '\xobs1'],'xobs1');
gif_maker_siman(R)
