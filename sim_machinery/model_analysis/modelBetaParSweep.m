%function [] = modelProbs(x,m,p,varo,R)
% R = simannealsetup_CSD_app()
clear ;close all
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\parsaves\prelim\naive.mat')
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\parsaves\prelim\modsetup_Rmu.mat')


R.IntP.dt = .0005;
R.obs.csd.df = 0.45;
R.obs.csd.reps = 12;

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = round(R.IntP.tend/R.IntP.dt);
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

% R.IntP.Utype = 'DCM_Str_Innov';
m.uset.p = p;
m.uset.p.covar = eye(m.m);
m.uset.p.scale = 1/1e-9; %.*R.IntP.dt;
u = innovate_timeseries(R,m);

% p = j;
x = m.x;

for jj = 1:1
%     ppm.increment();
    %% Simulate New Data
    % Integrate in time master fx function
    %         xsims = eval([R.IntP.intFx R.IntP.intFxArgs]);
    xsims = R.IntP.intFx(R,x,u,p,m);
    if isfield(R.obs,'obsFx') % Run Observer function
        xsims = R.obs.obsFx(xsims,m,p,R);
    end
figure; plot(xsims(1:6,:)'); shg
    if isfield(R.obs,'transFx') % Run Data Transform
        %% Construct CSD and compare to data
        %             fx = R.obs.transFx;
        %             [~,feat_sim] = fx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,10,R);
        [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,7,R);
    else
        feat_sim = xsims; % else take raw time series
    end
    % Now using NRMSE
    r2mean  = R.IntP.compFx(R,feat_sim);
    R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
%     r2rep{jj} = r2mean;
%     par_rep{jj} = pnew;
    disp(jj); % 
end

