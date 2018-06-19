function R = simannealsetup_NPD()
R.d = clock;

R.projectn = 'Rat_NPD';

if strmatch(getenv('computername'),'SFLAP-2')
    R.rootn = ['C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\sim_machinery';
else
    R.rootn = ['C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery';
end
addpath(genpath(R.rootn))
addpath(genpath(R.rootm))
R.filepathn = [R.rootn 'data\storage'];
R.data.datatype = 'CSD'; %%'NPD'

R.frqz = [8:.2:48];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MTX','STN','GPe','STR'};
R.chsim_name = {'MTX','STN','GPe','STR','GPi','THAL'};
R.out.tag = 'NPD_ABC';
% Set SimAn Parameters
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.A','.C'}; %,'.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.SimAn.pOptBound = [-5 3];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchN = 100;
R.SimAn.Tm = 0.8; % Initial temperature
R.SimAn.alpha = 0.98; % alpha increment
R.SimAn.rep = 56; % Repeats per temperature
R.SimAn.rtol_repeat = 0.85;
R.SimAn.rtol_converge = 0.95;
R.SimAn.ntol = 15;
R.SimAn.gradtol = [0.075 0.05];
R.SimAn.saveout = 'xobs1';
R.SimAn.maxdev = 3;
R.SimAn.jitter = 4;
R.SimAn.dSkew = 0.05;
R.SimAn.dPrec = 0.05;
R.SimAn.minRank = 38; %40;

% Set simulation parameters
% R.IntP.intFx = @stepintegrator_delay_Vec;
R.IntP.intFx = @spm_fx_compile;
R.IntP.intFxArgs = '(R,x,u,m,p)';
R.IntP.compFx= @compareData_100717;
R.IntP.compFxArgs = '(R,feat_sim)';

R.IntP.dt = .00085;
% R.IntP.tend = 10;
% R.IntP.nt = R.IntP.tend/R.IntP.dt;
% R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.n_Vet);
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

% Set Observer variables
R.obs.obsFx = @observe_data;
R.obs.obsFxArgs = '(xsims,m,pnew,R)';
R.obs.transFx = @constructCSDMat; %% @constructNPDMat;
R.obs.transFxArgs = '(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,[],R)';

R.obs.brn =1; % 2; % burn in time
R.obs.norm = 'False';
R.obs.csd.ztranscsd = 'True'; % z-transform CSDs
R.obs.csd.abovezero = 'True'; % Bring above zero

% desired freq res:
R.obs.csd.df = 0.45;
R.obs.csd.reps = 18;

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

% OR specify fft of 2^N:
% R.obs.csd.pow2_exp = 8;
% R.obs.csd.pow2_sim = 11;

R.obs.states = [7 9 11 13];
% LF = zeros(6);
% LF(1,1) = 0.1;
% LF(2,2) = 0.05;
% LF(3,3) = 0.05;
% LF(4,4) = 0.05;
% LF(5,5) = 0.05;
% LF(6,6) = 0.05;
%  LF = [1 0.05 0.05 0.05 0.05 0.05];
%         LF = [10 10 10 10 10 10];
LF = [20 2 2 2 2 2]; % Fit visually

R.obs.LF = LF;
R.obs.mixing = 0.1;

R.obs.gainmeth = {'leadfield','unitvar'};  %'mixing'
R.objfx.feattype = 'complex'; %%'ForRev'; %
R.objfx.specspec = 'cross'; %%'auto'; % which part of spectra to fit

R.plot.outFeatFx = @csdplotter_220517; %%@npdplotter_110717;
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;
R.plot.gif.delay = 0.3;
R.plot.gif.start_t = 1;
R.plot.gif.end_t = 1;
R.plot.gif.loops = 2;





