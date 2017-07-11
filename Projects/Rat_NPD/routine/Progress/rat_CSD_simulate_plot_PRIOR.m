%% RAT DATA- SIM ANNEAL PROJECT
% This script will fit BvW's BGM/MC model using simulated annealing, from
% the real part of the CSD for the group average OFF spectra
% TO DO:
% Extrinsic Inputs arent signed - all summed together. 
clear ; close all
R = simannealsetup;

%% Prepare the data

% prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\average_rat_lesion'])
% Set data as working version
ftdata = ftdata_lesion;
data = ftdata.trial{1}; fsamp = ftdata.fsample;
% Normalise Data
for i = 1:4
    xtmp = data(i,:);
    xtmp = (xtmp-mean(xtmp))/std(xtmp);
    data(i,:) = xtmp;
end
% compute CSD data features
[F_data,meancsd_data] = constructCSDMat(data,R.chloc_name,ftdata.label',fsamp,9,R.frqz);
% Plot CSD
csdplotter(meancsd_data,[],F_data,R)

% Set simulation parameters
dt = .001;
tend = 60;
nt = tend/dt;
tvec = linspace(0,tend,nt);

%% SIMULATION OF TIME SERIES
% 1 MMC
% 2 STR
% 3 GPE
% 4 STN
% 5 GPI
% 6 THAL

load('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\priors\rat_inversion\rat_inversion_params.mat')
clear u
[x,p,m] = setup_new_priors(P,M); % Assign Priors to correct structure
% setup exogenous noise
for src = 1:6
%     u{src} = rand(nt,m(1).nstates(src))/64;
    for j = 1:m(1).nstates(src)
        ux(:,j) = pinknoise(nt)/8;
    end
    u{src} = ux;
    clear ux
end

% Integrate in time master fx function
xstore = zeros(6,nt);
for t = 1:nt
    [x xobs] = simAn_master_fx_bgc(x,u,p,m,dt,t);
    xstore(:,t) = xobs';
    t
end

%% Process Simulated Data
dobs = [1 4 3 2]; % simulated data matching observed
% Normalise
for i = 1:length(dobs)
    xtmp = xstore(dobs(i),:);
    xtmp = (xtmp-mean(xtmp))/std(xtmp);
    xsims(i,:) = xtmp;
end
%% Construct CSD
[F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,ftdata.label',1/dt,10,R.frqz);
% Plot CSD
csdplotter([],meancsd_sim,F_sim,R)
