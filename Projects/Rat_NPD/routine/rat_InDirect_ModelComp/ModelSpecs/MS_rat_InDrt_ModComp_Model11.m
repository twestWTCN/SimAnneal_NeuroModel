function [R p m uc] = MS_rat_InDrt_ModComp_Model11(R)
% Model 11 - Hyperdirect/M2 Flat
m.m = 4; % # of sources
m.x = {[0 0 0 0 0 0 0 0] [0 0]  [0 0] [0 0]}; % Initial states
m.Gint = [14 1 1 1 1];
m.Tint = [4 1 1 1 1];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[0 0 0 0 0 0 1 0] [1 0] [1 0] [1 0]};
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
m.uset.p.scale = 5e1; %.*R.InstP.dt;
uc = innovate_timeseries(R,m);


%% Prepare Priors
% 1 MMC
% 2 STR
% 3 GPE
% 4 STN

% Excitatory connections
p.A{1} =  repmat(-32,m.m,m.m);
p.A{1}(2,1) = 0; % MMC -> STR
p.A{1}(3,4) = 0; % STN -> GPe
p.A{1}(4,1) = 0; % MMC -> STN (hyperdirect)
p.A_s{1} = repmat(1,m.m,m.m);

p.A{2} =  repmat(-32,m.m,m.m);
p.A{2}(3,2) = 0; % STR -| GPe
p.A{2}(4,3) = 0; % GPe -| STN
p.A_s{2} = repmat(1,m.m,m.m);

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.5,size(p.C));

% Leadfield
p.obs.LF = [1 1 1 1];
p.obs.LF_s = repmat(0.2,size(p.obs.LF));

p.obs.mixing = [1]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0,size(p.obs.mixing));

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.25,size(p.D));

% Sigmoid transfer for connections
p.S = [0 0];
p.S_s = [0.2 0.2];

% time constants and gains
for i = 1:m.m
    if i == 1
        prec = 1;
    else
        prec = 1;
    end
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(prec,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(prec,size(p.int{i}.G));
    p.int{i}.S = zeros(1);
    p.int{i}.S_s = repmat(prec,size(p.int{i}.S));
    %     p.int{i}.BT = zeros(1,m.Tint(i));
    %     p.int{i}.BT_s = repmat(prec,size(p.int{i}.T));
end

% MMC Auto Flat Priors
prec= 0.25;
p.int{1}.T = [-0.14 0.15 -0.74 0.11];
p.int{1}.T_s = repmat(prec,size(p.int{1}.T));
p.int{1}.G =[0.05 0.31 0.24 0.08 0.66 0.08 -0.15 -0.70 0.01 -0.12 0.60 0.37 0.72 0.09];
p.int{1}.G_s = repmat(prec,size(p.int{1}.G));
p.int{1}.S = 0.98;
p.int{1}.S_s =  repmat(prec,size(p.int{1}.S));

% % %% SPECIFIC TO STN/GPe Resonator
% % % GPe
% % p.int{3}.T = 0.67;
% % p.int{3}.T_s = prec;
% % p.int{3}.G = -0.64;
% % p.int{3}.G_s = prec;
% % p.int{3}.S = 0.16;
% % p.int{3}.S_s = prec;
% % % STN
% % p.int{4}.T = 1.5;
% % p.int{4}.T_s = prec;
% % p.int{4}.G = 2.25;
% % p.int{4}.G_s = prec;
% % p.int{4}.S = 0.09;
% % p.int{4}.S_s = prec;