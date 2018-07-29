%% RAT DATA- SIM ANNEAL PROJECT
%%%%%%%%%%%%%%%%%%
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% 1) Take all sims that are above criteria and include!
% 2) Set Tm and anneal to be optimal - adjust to the gradient of the decent
%%%%%%%%%%%%%%%%%%%%%%%%
% simAnnealAddPaths()
clear ; close all
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery'))
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD'))
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools\')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot\')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
% addpath('C:\spm12')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
%
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
rng(5324)

%% Set Parameters of the Routine
R = simannealsetup_NPD_MMC;
R.d = clock; % clock

%% Prepare the data
% prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\NPD_paper_RatNPD_150618.mat']);
NPDmat = fyA;
load([R.filepathn '\nsPow_paper_RatNPD_150618.mat']);
nsPowmat = fyA;
load([R.filepathn '\frq_paper_RatNPD_150618.mat']);
F_data = fxA(:,1);
meannpd_data = [];
% condsel = [1 2];
R.condnames = {'ON'};
R.Bcond = 1;
condsel = 1;
chsel = [1];
for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                Ftmp = F_data;
                Pxy = abs(log10(mean(vertcat(nsPowmat{chsel(i),condsel(C),:}),1)));
                Pxy(Ftmp>48 & Ftmp<52) = [];
                Ftmp(Ftmp>48 & Ftmp<52) = [];
                [xCalc yCalc b Rsq] = linregress(log10(Ftmp),Pxy');
                Pxy = Pxy-yCalc';
                Pxy = interp1(Ftmp,Pxy,R.frqz);
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy - min(Pxy);
                Pxy = Pxy.* gausswin(length(Pxy),1)';% tukeywin(length(Pxy),0.1)';
                 plot(R.frqz,Pxy)
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                for k = 1:size(NPDmat,3)
                    Ftmp = F_data;
                    Cxy = mean(horzcat(NPDmat{chsel(i),chsel(j),k,condsel(C),:})',1);
                    Cxy = Cxy*2;
                    Cxy(Ftmp>48 & Ftmp<52) = [];
                    Ftmp(Ftmp>48 & Ftmp<52) = [];
                    Cxy = interp1(Ftmp,Cxy,R.frqz);
                    Cxy = Cxy.*tukeywin(length(Cxy),0.1)'; %%NPD_sim_n(i,j,1,:)
%                     plot(R.frqz,Cxy)
                    meannpd_data(C,i,j,k,:) = Cxy;
                end
            end
        end
    end
end
% Set data as working version
R.data.feat_emp = meannpd_data;
% squeeze(meannpd_data(1,1,1,1,:))
R.data.feat_xscale = R.frqz;

% Plot CSD
if strcmp('CSD',R.data.datatype)
    csdplotter_220517({meannpd_data},[],F_data,R)
elseif strcmp('NPD',R.data.datatype)
    npdplotter_110717({meannpd_data},[],R.frqz,R,[],[])
end

%% Prepare Model
m.m = 1; % # of sources
m.x = {[0 0 0 0 0 0 0 0]}; % Initial states
m.Gint = [14];
m.Tint = [4];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[0 0 0 0 0 0 1 0]};
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
% 5 GPI
% 6 THAL

% Excitatory connections
p.A{1} = -32;
p.A_s{1} = 0;
p.A{2} = -32;
p.A_s{2} = 0;
% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.5,size(p.C));

% Leadfield
p.obs.LF = [1 1];
p.obs.LF_s = repmat(0.2,size(p.obs.LF));

p.obs.mixing = [1]; %zeros(size(R.obs.mixing));
p.obs.mixing_s = repmat(0,size(p.obs.mixing));

% Delays
p.D = -32;
p.D_s = 0;

% Sigmoid transfer for connections
p.S = [0 0];
p.S_s = [0.2 0.2];

% time constants and gains
for i = 1:m.m
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(1,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(1,size(p.int{i}.G));
    p.int{i}.S = zeros(1);
    p.int{i}.S_s = repmat(1,size(p.int{i}.S));
    p.int{i}.BT = zeros(1,m.Tint(i));
    p.int{i}.BT_s = repmat(1,size(p.int{i}.T));   
end
pold = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R.out.dag = 'NPD_MMC_ON'; % 'All Cross'
R.out.tag = 'NPD_MMC_ON';
% R.out.dag = '180718_COND1_auto'; % cross only
% R.plot.save = 'True';
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
% m = varo;
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
% xobs1 = varo;
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
% p = xobs1.Mfit.Pfit;
% p.A = pold.A;
% p.A{1}(2,4) = 0; % STN to STR
% p.D(2,4) = 0;
% A_s = repmat(0,size(A));
% p.A_s{1} = A_s;
% p.A_s{2} = A_s;
% 
% % parBank = varo;
R.SimAn.rep =40; %512; %96; %512; % Repeats per temperature
R.SimAn.Tm = 1;
R.SimAn.jitter = 1;
R.SimAn.searchN = 200;
R = setSimTime(R,12);
R.objfx.specspec = 'auto';
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C'}; %,'.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.Bcond = 0;
parBank = [];
R.SimAn.copout = [4 6];
[xcross_OFF] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);

