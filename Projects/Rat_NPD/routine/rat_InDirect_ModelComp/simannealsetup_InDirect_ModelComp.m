function R = simannealsetup_InDirect_ModelComp(R)
% 
R.projectn = 'Rat_NPD';
R.out.tag = 'InDrt_ModCompRev2';
R.out.oldtag = 'InDrt_ModCompRev2';

switch getenv('computername')
    case 'SFLAP-2'
    R.rootn = ['C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\sim_machinery';
    R.BBA_path = 'C:\Users\Tim\Documents\Work\GIT\Cortical_Parkinsons_Networks\BetaBursts';
    case 'FREE'
    R.rootn = ['C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery';
    R.BBA_path = 'C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks\BetaBursts';
    case 'DESKTOP-94CEG1L'
    R.rootn = ['C:\Users\timot\Documents\GitHub\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\timot\Documents\Github\SimAnneal_NeuroModel\sim_machinery';
    R.BBA_path = 'C:\Users\timot\Documents\GitHub\BurstToolbox';
end
% addpath(genpath(R.rootn))
% addpath(genpath(R.rootm))
% 
%% DATA SPECIFICATION
R.filepathn = [R.rootn 'data\storage'];
R.data.datatype = 'NPD'; %%'NPD'
R.frqz = [6:.2:68];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MMC','STR','GPE','STN'};
R.chsim_name = {'MMC','STR','GPE','STN'};
R.condnames = {'OFF'};
% Spectral characteristics
R.obs.csd.df = 0.5;
R.obs.csd.reps = 32; %36; %96;

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @spm_fx_compile_120319;
R.IntP.compFx= @compareData_100717;

R.IntP.dt = .001;
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

%% OBSERVATION 
% observation function
R.obs.obsFx = @observe_data;
R.obs.gainmeth = {'unitvar','boring'}; %,'submixing'}; %,'lowpass'}; ,'leadfield' %unitvar'mixing'
R.obs.glist =0; %linspace(-5,5,12);  % gain sweep optimization range [min max listn] (log scaling)
R.obs.brn =3; % 2; % burn in time
LF = [1 1 1 1]*10; % Fit visually and for normalised data
R.obs.LF = LF;
% % (precompute filter)
% % fsamp = 1/R.IntP.dt;
% % nyq = fsamp/2;
% % Wn = R.obs.lowpass.order/nyq;
% % R.obs.lowpass.fwts = fir1(R.obs.lowpass.order,Wn);

% Data Features
% fx to construct data features
R.obs.transFx = @constructNPDMat_190618;
% These are options for transformation (NPD)
R.obs.trans.logdetrend =0;
R.obs.trans.norm = 1;
R.obs.trans.gauss = 0;
R.obs.logscale = 0;


%% OBJECTIVE FUNCTION
R.objfx.feattype = 'ForRev'; %%'ForRev'; %
R.objfx.specspec = 'cross'; %%'auto'; % which part of spectra to fit

%% OPTIMISATION
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.C','.A','.D'}; %,'.int{src}.S','.A',,'.int{src}.BG','.int{src}.S','.S','.D','.obs.LF'};  %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',
R.SimAn.pOptBound = [-12 12];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);
R.SimAn.searchMax = 200;
R.SimAn.convIt = 1e-4;
R.SimAn.rep = 720; %128; %512; % Repeats per temperature
R.SimAn.saveout = 'xobs1';
R.SimAn.jitter = 1; % Global precision
%% PLOTTING
R.plot.outFeatFx = @npdplotter_110717; %%@;csdplotter_220517
R.plot.save = 'False';
R.plot.distchangeFunc = @plotDistChange_KS;
R.plot.gif.delay = 0.3;
R.plot.gif.start_t = 1;
R.plot.gif.end_t = 1;
R.plot.gif.loops = 2;




%% OLD PARAMETERS
% R.obs.transFxArgs = '(xsims_gl{gl},R.chloc_name,R.chsim_name,1/R.IntP.dt,R.obs.SimOrd,R)';
% % R.obs.norm = 'False';
% % R.obs.csd.ztranscsd = 'False'; % z-transform CSDs
% % R.obs.csd.abovezero = 'False'; % Bring above zero
% % R.SimAn.starttemp = 2;
% % R.obs.mixing = [0.005 0.05];
% % R.obs.lowpass.cutoff = 80;
% % R.obs.lowpass.order = 80;
% % R.SimAn.maxdev = 12;
% % R.SimAn.starttemp = 2;
% % R.SimAn.alpha = 0.98; % alpha increment


