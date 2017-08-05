%function [] = modelProbs(x,m,p,varo,R)
% R = simannealsetup_CSD_app()
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\parsaves\prelim\sim_ABC_output_010817_d.mat')
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\parsaves\prelim\modsetup_Rmu.mat')


R.IntP.dt = .0004;
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

u = innovate_timeseries(R,m);

p = j;


for jj = 1:1
%     ppm.increment();
    %% Simulate New Data
    % Integrate in time master fx function
    %         xsims = eval([R.IntP.intFx R.IntP.intFxArgs]);
    xsims = R.IntP.intFx(R,x,u,p,m);
    
    if isfield(R.obs,'obsFx') % Run Observer function
        xsims = R.obs.obsFx(xsims,m,p,R);
    end
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

