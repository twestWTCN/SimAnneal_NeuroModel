%% RAT DATA- SIM ANNEAL PROJECT
%%%%%%%%%%%%%%%%%%
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
%
%%%%%%%%%%%%%%%%%%%%%%%%

clear ; close all
addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery'))
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools\')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot')
addpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\model_fx')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
addpath('C:\spm12')
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\ParforProgMon'))
addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\boundedline-pkg'))

rng(623)

%% Set Parameters of the Routine
R = simannealsetup_comparNPD_Granger_170118;
R.d = clock; % clock

load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\NPD_Granger_compar\Data\3Node_MVAR_Data.mat')
R.obs.obsstates = [1 2 3];
R.obs.logdetrend = 0;
[F meannpd_data] = constructNPDMat_190618(data.trial,R.chloc_name,R.chloc_name,data.fsample/2,9,R)

for z= 1:3
    for i = 1:3
        for j = 1:3
            if i~=j
                meannpd_data(1,z,i,j,:) = 0.75.*meannpd_data(1,z,i,j,:);
            end
        end
    end
end
R.data.feat_emp = meannpd_data;
R.data.feat_xscale = R.frqz;

%% Prepare Model
m.m = 3; % # of sources
m.x = {[0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0] [0 0 0 0 0 0 0 0]}; % Initial states
m.Gint = [8 8 8] ; % Priors
m.Tint = [4 4 4];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
m.dipfit.model(1).source = 'MMC';
m.dipfit.model(2).source = 'MMC';
m.dipfit.model(3).source = 'MMC';


m.outstates = {[0 0 0 0 0 0 1 0] [0 0 0 0 0 0 1 0] [0 0 0 0 0 0 1 0] };
R.obs.outstates = find([m.outstates{:}]);
for i=1:numel(R.chloc_name)
    R.obs.obsstates(i) = find(strcmp(R.chloc_name{i},R.chsim_name));
end

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
% m.uset.p.scale = 5e-2; %.*R.InstP.dt;
rng(8453)
% uc = innovate_timeseries(R,m);

%% Prepare Priors
% 1 X1
% 2 X2
% Excitatory connections
A = repmat(-32,m.m,m.m);
% A(2,1) = -12;
% A(2,1) = 1;
% A(3,1) = 1;

p.A{1} = A;
A_s = repmat(1,size(A));
p.A_s{1} = A_s;

% Inhbitory connections
A = repmat(-32,m.m,m.m);
% A(2,1) = 6;
p.A{2} = A;
A_s = repmat(1,size(A));
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
p.D = repmat(0,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.1,size(p.D));

% Sigmoid transfer for connections
p.S = [0 0];
p.S_s = [0.1 0.1];

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
R.out.tag = 'comparNPD_Granger';
R.out.dag = '100718';

% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat']);
% m = varo;
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat']);
% xobs1 = varo;
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat']);
% parBank = varo;
parBank = [];

R.SimAn.jitter = 1;
m.uset.p.scale = 1e-4;
R.SimAn.searchN = 50;


% Now do Cross
R.plot.save = 'True';
A = repmat(1,m.m,m.m);
% A(2,1) = 0;
% A(3,1) = 0;
% A(1,1) = 0;
% A(1,3) = 0;

p.A{1} = A;
A_s = repmat(0.5,size(A));
p.A_s{1} = A_s;

A = repmat(-32,m.m,m.m);
p.A{2} = A;
A_s = repmat(0.5,size(A));
p.A_s{2} = A_s;
%
parBank = [];
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.S','.C','.A'}; %,'.S','.obs.LF'}; ,'.int{src}.S','.S' %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',

R.out.dag = [R.out.dag 'cross'];
R.objfx.specspec = 'cross';
R.SimAn.searchN = 100;
R = setSimTime(R,24);

for i = 1:6
    %     R.out.tag = [R.out.tag num2str(i)];
    if i>1
        p = xobs1.Mfit.Pfit;
        R.plot.save = 'False';
        %         R.SimAn.searchN = R.SimAn.searchN * 0.5;
        R.SimAn.Tm = 0.75; % R.SimAn.Tm -0.1;
    end
    [xobs1 parBank] = SimAn_ABC_110817(m.x,[],p,m,R,parBank);
end
% folname = ['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\parfits\' sprintf('%d',[d(1:3)])];

gif_maker_siman(R)
