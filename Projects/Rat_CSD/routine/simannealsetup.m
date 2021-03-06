function R = simannealsetup()
R.d = clock;

R.rootn = 'C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\';
R.projectn = 'Rat_CSD';
R.filepathn = [R.rootn R.projectn '\data\storage'];
R.frqz = [8:.2:40];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MTX','STN','GPe','STR'};
R.chsim_name = {'MTX','STN','GPe','STR','GPi','THAL'};
R.out.tag = 'autost';
% Set SimAn Parameters 
R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.C','.D','.A','.obs.LF','.obs.mixing'};
R.SimAn.pOptBound = [-6 6];

R.SimAn.searchN = 50;
R.SimAn.Tm = 1; % Initial temperature
R.SimAn.alpha = 0.98; % alpha increment
R.SimAn.rep = 16; % Repeats per temperature
R.SimAn.rtol_repeat = 0.85;
R.SimAn.rtol_converge = 0.95;
R.SimAn.ntol = 15;
R.SimAn.gradtol = [0.075 0.05];
R.SimAn.saveout = 'xobs1';
R.SimAn.maxdev = 6;
R.SimAn.jitter = 4;
R.SimAn.dSkew = 30;
R.SimAn.dPrec = 25;

R.SimAn.opPar = {'T','G','C','D','A','LF','mix'};
% Set simulation parameters
R.IntP.intFunc = @stepintegrator_delay;
R.IntP.dt = .0001;
R.IntP.tend = 8;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

% Set Observer variables
R.obs.fx = @observe_data;
R.obs.brn =1; % 2; % burn in time
R.obs.norm = 'False';
R.obs.csd.ztranscsd = 'True'; % z-transform CSDs
R.obs.csd.abovezero = 'True'; % Bring above zero

% desired freq res:
R.obs.csd.df = 1.5;
R.obs.csd.reps = 22;

fsamp = 1/R.IntP.dt;    % sample rate
N = floor(fsamp/R.obs.csd.df);  % segment length
R.IntP.tend = (N*R.obs.csd.reps)/fsamp; % Target time
disp(sprintf('The simulation length is %.2f seconds',R.IntP.tend))
disp(sprintf('The simulation df is %.2f Hz',R.obs.csd.df))

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
 LF = [8 4 1.2 20 10 10]; % Fit visually

R.obs.LF = LF;
R.obs.mixing = 0.1;

R.obs.gainmeth = {'unitvar','mixing'};  %'leadfield'
R.objfx.feattype = 'complex';
R.objfx.specspec = 'cross'; % which part of spectra to fit

R.plot.gif.delay = 0.3;
R.plot.gif.start_t = 1;
R.plot.gif.end_t = 1;
R.plot.gif.loops = 2;

if strmatch(getenv('COMPUTERNAME'),'SFLAP-2')
    R.rootdir = 'C:\Users\Tim\Documents\Work\LitvakProject\SimAnneal_NeuroModel\';
elseif strmatch(getenv('COMPUTERNAME'),'FREE')
    R.rootdir = 'C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\';
end



