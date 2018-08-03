% THIS IS THE STN/GPE with Independent M2

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
% %
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
rng(23123)

%% Set Parameters of the Routine
R = simannealsetup_STN_GPe_M2_ModelComp;
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
chsel = [1 2 4]';

for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                Ftmp = F_data;
                Pxy = abs(log10(mean(vertcat(nsPowmat{chsel(i),condsel(C),:}),1)));
                %                 Pxy = Pxy(Ftmp>=R.frqz(1) & Ftmp<=R.frqz(end));
                %                 Ftmp = Ftmp(Ftmp>=R.frqz(1) & Ftmp<=R.frqz(end));
                Pxy(Ftmp>48 & Ftmp<52) = [];
                Ftmp(Ftmp>48 & Ftmp<52) = [];
                [xCalc yCalc b Rsq] = linregress(log10(Ftmp),Pxy');
                Pxy = Pxy-yCalc';
                Pxy = interp1(Ftmp,Pxy,R.frqz);
                %                     plot(R.frqz,Pxy); hold on
                f = fit(R.frqz',Pxy','gauss3');
                Pxy = f(R.frqz)';
                Pxy = 10.^Pxy;
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy - min(Pxy);
                %                 Pxy = Pxy.*tukeywin(length(Pxy),0.6)';
                %                     plot(R.frqz,Pxy); hold on
                %                     plot(R.frqz,10.^(Pxy));
                meannpd_data(C,i,j,1,:) = Pxy;
            else
                if i == 1 || j == 1
                    for k = 1:size(NPDmat,3)
                        Cxy = exp(-4).*randn(size(R.frqz));
                        meannpd_data(C,i,j,k,:) = Cxy;
                    end
                else
                    for k = 1:size(NPDmat,3)
                        Ftmp = F_data;
                        Cxy = mean(horzcat(NPDmat{chsel(i),chsel(j),k,condsel(C),:})',1);
                        Cxy = Cxy*1.5;
                        Cxy(Ftmp>48 & Ftmp<52) = [];
                        Ftmp(Ftmp>48 & Ftmp<52) = [];
                        Cxy = interp1(Ftmp,Cxy,R.frqz);
                        %                     Cxy = Cxy.*tukeywin(length(Cxy),0.6)'; %%NPD_sim_n(i,j,1,:)
                        %                     plot(R.frqz,Cxy); hold on
                        f = fit(R.frqz',Cxy','gauss3');
                        Cxy = f(R.frqz)';
                        %                     plot(R.frqz,Cxy)
                        meannpd_data(C,i,j,k,:) = Cxy;
                        %                     close all
                    end
                end
            end
        end
    end
end

%% Prepare the data
% prepareratdata_group(R.rootn,R.projectn);
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
m.m = 3; % # of sources
m.x = {[0 0 0 0 0 0 0 0] [0 0]  [0 0]}; % Initial states
m.Gint = [14 1 1];
m.Tint = [4 1 1];
m.n =  size([m.x{:}],2); % Number of states
% These outline the models to be used in compile function
for i = 1:numel(R.chsim_name)
    m.dipfit.model(i).source = R.chsim_name{i};
end

m.outstates = {[0 0 0 0 0 0 1 0] [1 0] [1 0]};
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
% 3 GPE
% 4 STN

% Excitatory connections
p.A{1} =  repmat(-32,m.m,m.m);
p.A{1}(2,3) = 0;
p.A_s{1} = repmat(1,m.m,m.m);
p.A_s{1}(2,3) = 1.5;

p.A{2} =  repmat(-32,m.m,m.m);
p.A{2}(3,2) = 0;
p.A_s{2} = repmat(1,m.m,m.m);
p.A_s{2}(3,2) = 1.5;

% Connection strengths
p.C = zeros(m.m,1);
p.C_s = repmat(0.5,size(p.C));

% Leadfield
p.obs.LF = [1 1];
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
        prec = 0.5;
    else
        prec = 1.5;
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
pold = p;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R.out.dag = 'NPD_ModComp'; % 'All Cross'
R.out.tag = 'Mod1';
R.obs.trans.norm = 1;
R.obs.logdetrend =1;

R.SimAn.rep =64; %512; %96; %512; % Repeats per temperature
R.SimAn.Tm = 1;
R.SimAn.jitter = 1.5;
R.SimAn.searchN = 200;
R = setSimTime(R,18);
R.objfx.specspec = 'cross';
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A',}; %,'.D','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.Bcond = 0;
parBank = [];
R.SimAn.copout = [2 3];
[xcross_OFF] = SimAn_ABC_110817(m.x,uc,p,m,R,parBank);

