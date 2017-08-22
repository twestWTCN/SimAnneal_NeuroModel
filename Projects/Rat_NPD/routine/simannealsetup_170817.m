function R = simannealsetup_170817()
% 
R.projectn = 'Rat_NPD';
R.out.tag = 'CSD_neatmodel';

if strmatch(getenv('computername'),'SFLAP-2')
    R.rootn = ['C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\sim_machinery';
else
    R.rootn = ['C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery';
end
addpath(genpath(R.rootn))
addpath(genpath(R.rootm))


%% DATA SPECIFICATION
R.filepathn = [R.rootn 'data\storage'];
R.data.datatype = 'CSD'; %%'NPD'
R.frqz = [6:.2:85];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MTX','STN','GPe','STR'};
R.chsim_name = {'MTX','STN','GPe','STR','GPi','THAL'};

%% OBSERVATION 
R.obs.obsFx = @observe_data;
R.obs.obsFxArgs = '(xsims,m,pnew,R)';
R.obs.gainmeth = {'unitvar','leadfield','submixing'}; %,'lowpass'};  %unitvar'mixing'

R.obs.transFx = @constructCSDMat; %% @constructNPDMat;
R.obs.transFxArgs = '(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,[],R)';
R.obs.brn =1; % 2; % burn in time
R.obs.norm = 'False';
R.obs.csd.ztranscsd = 'False'; % z-transform CSDs
R.obs.csd.abovezero = 'False'; % Bring above zero
% desired freq res:
R.obs.csd.df = 0.45;
R.obs.csd.reps = 12;
R.obs.states = [7 9 11 13];
% LF = [1 0.008 0.005 0.005 0.005 0.005]; % for non-normalised
% LF = [1 1 1 1 1 1]; % Fit visually and for normalised data
  LF = [0.8 0.5 0.5 0.5 0.5 0.5]; % Fit visually and for normalised data
R.obs.LF = LF;
R.obs.mixing = [0.005 0.05];
R.obs.lowpass.cutoff = 80;
R.obs.lowpass.order = 80;

%% INTEGRATION
R.IntP.intFx = @spm_fx_compile_180817;
R.IntP.intFxArgs = '(R,x,u,m,p)';
R.IntP.compFx= @compareData_100717;
R.IntP.compFxArgs = '(R,feat_sim)';

R.IntP.dt = .0005;
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

% (precompute filter)
fsamp = 1/R.IntP.dt;
nyq = fsamp/2;
Wn = R.obs.lowpass.order/nyq;
R.obs.lowpass.fwts = fir1(R.obs.lowpass.order,Wn);


%% OBJECTIVE FUNCTION
R.objfx.feattype = 'complex'; %%'ForRev'; %
R.objfx.specspec = 'auto'; %%'auto'; % which part of spectra to fit

%% OPTIMISATION
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.S','.A','.D','.C','.obs.mixing','.obs.LF'}; %,'.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchN = 100;
R.SimAn.Tm = 0.8; % Initial temperature
R.SimAn.alpha = 0.97; % alpha increment
R.SimAn.rep = 128; % Repeats per temperature
R.SimAn.saveout = 'xobs1';
% R.SimAn.maxdev = 12;
R.SimAn.jitter = 1;

%% PLOTTING
R.plot.outFeatFx = @csdplotter_220517; %%@npdplotter_110717;
R.plot.save = 'True';
R.plot.distchangeFunc = @plotDistChange_KS;
R.plot.gif.delay = 0.3;
R.plot.gif.start_t = 1;
R.plot.gif.end_t = 1;
R.plot.gif.loops = 2;





