%% Simulation analyses
clear all
% Run optimized parameters
R = simannealsetup;
load([R.rootdir 'Projects\Rat_CSD\sim_analyses\parsaves\xobs1.mat'])
load([R.rootdir 'Projects\Rat_CSD\sim_analyses\parsaves\fpar_050717.mat'])
xstore = stepintegrator_delay(R,x,u,m,xobs1.out.P);
% Run Observer function
xsims = observe_data(xstore,m,xobs1.out.P,R);
%% Construct CSD and compare to data
[F_sim,meancsd_sim] = constructCSDMat(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,[],R);
% Now using NRMSE
r2mean = compareCSD_240517(R.data.meancsd_data,meancsd_sim,R);

csdplotter_220517({R.data.meancsd_data},{meancsd_sim},R.frqz,R)