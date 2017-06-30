function R = simannealsetup()
R.rootn = 'C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\';
R.projectn = 'Rat_CSD';
R.filepathn = [R.rootn R.projectn '\data\storage'];
R.frqz = [10:.2:80];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MTX','STN','GPe','STR'};
R.chsim_name = {'MTX','STN','GPe','STR','GPi','THAL'};

% Set SimAn Parameters 
R.SimAn.searchN = 100;
R.SimAn.Tm = 1; % Initial temperature
R.SimAn.alpha = 0.99; % alpha increment
R.SimAn.rep = [16 8]; %24; % Repeats per temperature
R.SimAn.rtol_repeat = 0.85;
R.SimAn.rtol_converge = 0.95;
% Set simulation parameters
R.IntP.dt = .001;
R.IntP.tend = 10;
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);
R.IntP.Utype = 'white_covar'; %'white_covar'; % DCM_Str_Innov

% Set Observer variables
R.obs.brn =1; % 2; % burn in time
R.obs.norm = 'False';
R.obs.csd.ztranscsd = 'True'; % z-transform CSDs
R.obs.csd.abovezero = 'True'; % Bring above zero
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
R.objfx.specspec = 'auto'; % which part of spectra to fit






