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
rng(421)

%% Set Parameters of the Routine
R = simannealsetup_NPD_Can_060718;
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
R.condnames = {'OFF'};
R.Bcond = 1;
condsel = 2;
for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(NPDmat,1)
        for j = 1:size(NPDmat,2)
            if i==j
                Ftmp = F_data;
                Pxy = abs(log10(mean(vertcat(nsPowmat{i,condsel(C),:}),1)));
                Pxy(Ftmp>48 & Ftmp<52) = [];
                Ftmp(Ftmp>48 & Ftmp<52) = [];
                [xCalc yCalc b Rsq] = linregress(log10(Ftmp),Pxy');
                Pxy = Pxy-yCalc';
                Pxy = interp1(Ftmp,Pxy,R.frqz);
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy - min(Pxy);
                Pxy = Pxy; %.*tukeywin(length(Pxy),0.25)';
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                for k = 1:size(NPDmat,3)
                    Ftmp = F_data;
                    Cxy = mean(horzcat(NPDmat{i,j,k,condsel(C),:})',1);
                    Cxy = Cxy*1.5;
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
m.m = 6; % # of sources
m.x = {[0 0 0 0 0 0 0 0]  [0 0]  [0 0]  [0 0]  [0 0]  [0 0]}; % Initial states
m.Gint = [14 1 1 1 1 1];
m.Tint = [4 1 1 1 1 1];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[0 0 0 0 0 0 1 0]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]};
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
A = repmat(-32,m.m,m.m);
A(2,1) = 0; % M1 to STR
A(4,1) = 0; % M1 to STN
A(3,4) = 0; % STN to GPe
A(5,4) = 0; % STN to GPi
A(1,6) = 0; % THAL to M1

p.A{1} = A;
A_s = repmat(1,size(A));
p.A_s{1} = A_s;

% Inhbitory connections
A = repmat(-32,m.m,m.m);
A(3,2) = 0; % STR to GPe
A(5,2) = 0; % STR to GPi
A(4,3) = 0; % GPe to STN
A(5,3) = 0; % GPe to GPi
A(2,3) = 0; % GPe to STR
A(6,5) = 0; % GPi to THAL
p.A{2} = A;
A_s = repmat(2,size(A));
p.A_s{2} = A_s;

% % p.B{1} = p.A{1}; p.B{2} = p.A{2};
% % p.B{1}(p.A{1}>-30) = 0; %zero!!
% % p.B{2}(p.A{2}>-30) = 0; %zero!!
% % p.B_s{1} = repmat(1,size(A));
% % p.B_s{2} = repmat(1,size(A));

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.2,size(p.C));

% Leadfield
p.obs.LF = zeros(size(R.obs.LF)).*0.8;
p.obs.LF_s = repmat(0.2,size(p.obs.LF));

% Delays
p.D = repmat(-32,size(p.A{1})).*~((p.A{1}>-32) | (p.A{2}>-32)) ;
p.D_s = repmat(0.2,size(p.D));

% Sigmoid transfer for connections
p.S = [0 0];
p.S_s = [0.2 0.2];

% time constants and gains
for i = 1:m.m
    p.int{i}.T = zeros(1,m.Tint(i));
    p.int{i}.T_s = repmat(0.5,size(p.int{i}.T));
    p.int{i}.G = zeros(1,m.Gint(i));
    p.int{i}.G_s = repmat(0.5,size(p.int{i}.G));
    p.int{i}.S = zeros(1);
    p.int{i}.S_s = repmat(0.5,size(p.int{i}.S));
    
    p.int{i}.BT = zeros(1,m.Tint(i));
    p.int{i}.BT_s = repmat(1.5,size(p.int{i}.T));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BANK
%  '180718_COND1_reauto' - Really nice auto fit!
%  '180718_COND1_full' - partially fitted to all(~40% R2)

 % Now fit Auto  
 R.out.dag = '180718_COND1_full'; % cross only
R.plot.save = 'True';
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
xobs1 = varo;
p = xobs1.Mfit.Pfit;
% load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
% parBank = varo;
parBank = [];
p.A{1}(3,4) = 0;
p.A{2}(4,3) = 0.5;
p.A_s{1}(3,4) = 2;
p.A_s{2}(4,3) = 2;

R.out.dag = '180718_COND1_full_run2'; % cross only
R.objfx.specspec = 'cross';
R = setSimTime(R,32);
R.SimAn.rep = 128;
R.SimAn.Tm = 1;
R.SimAn.jitter =1;
R.Bcond = 0;
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.A','.S','.C'}; %
[xcross_OFF] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);
 
    
% % %     
% % %     
% % % %% Now Find  ON
% % % R.condnames = {'ON'};
% % % condsel = 1;
% % % for C =1:numel(R.condnames)
% % %     %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
% % %     for i = 1:size(NPDmat,1)
% % %         for j = 1:size(NPDmat,2)
% % %             if i==j
% % %                 Ftmp = F_data;
% % %                 Pxy = abs(log10(mean(vertcat(nsPowmat{i,condsel(C),:}),1)));
% % %                 Pxy(Ftmp>48 & Ftmp<52) = [];
% % %                 Ftmp(Ftmp>48 & Ftmp<52) = [];
% % %                 [xCalc yCalc b Rsq] = linregress(log10(Ftmp),Pxy');
% % %                 Pxy = Pxy-yCalc';
% % %                 Pxy = interp1(Ftmp,Pxy,R.frqz);
% % %                 Pxy = (Pxy-mean(Pxy))./std(Pxy);
% % %                 Pxy = Pxy - min(Pxy);
% % %                 Pxy = Pxy; %.*tukeywin(length(Pxy),0.25)';
% % %                 meannpd_data(C,i,j,1,:) = Pxy;
% % %             else
% % %                 for k = 1:size(NPDmat,3)
% % %                     Ftmp = F_data;
% % %                     Cxy = mean(horzcat(NPDmat{i,j,k,condsel(C),:})',1);
% % %                     Cxy = Cxy*2;
% % %                     Cxy(Ftmp>48 & Ftmp<52) = [];
% % %                     Ftmp(Ftmp>48 & Ftmp<52) = [];
% % %                     Cxy = interp1(Ftmp,Cxy,R.frqz);
% % %                     Cxy = Cxy.*tukeywin(length(Cxy),0.1)'; %%NPD_sim_n(i,j,1,:)
% % %                     plot(R.frqz,Cxy)
% % %                     meannpd_data(C,i,j,k,:) = Cxy;
% % %                 end
% % %             end
% % %         end
% % %     end
% % % end
% % % % Set data as working version
% % % R.out.dag = '180718_COND2_auto'; % cross only
% % % R.data.feat_emp = meannpd_data;
% % % % % p = xcross_OFF.Mfit.Pfit;
% % %  for i = 1:m.m
% % %     p.int{i}.T_s = repmat(0.1,size(p.int{i}.T));
% % %     p.int{i}.G_s = repmat(0.1,size(p.int{i}.G));
% % %     p.int{i}.S_s = repmat(0.1,size(p.int{i}.S));
% % %     p.int{i}.BT_s = repmat(1,size(p.int{i}.T));
% % %  end
% % % R.Bcond = 1; % The condition which is B additive
% % % R.objfx.specspec = 'auto';
% % % R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.BT','.S','.C'}; %,'.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
% % % parBank = [];
% % % 
% % % [xauto_on] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);
% % % 
% % % 
% % % % folname = ['C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\outputs\parfits\' sprintf('%d',[d(1:3)])];
% % % % mkdir(folname)
% % % % save([folname '\xobs1'],'xobs1');
% % % gif_maker_siman(R)
