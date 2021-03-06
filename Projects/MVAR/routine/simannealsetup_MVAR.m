function R = simannealsetup_MVAR()
R.d = clock;
R.projectn = 'MVAR';
if strmatch(getenv('computername'),'SFLAP-2') 
    R.rootn = ['C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\sim_machinery';
else
    R.rootn = ['C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\Projects\' R.projectn '\'];
    R.rootm = 'C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery';
end
addpath(genpath(R.rootn))
addpath(genpath(R.rootm))


R.filepathn = [R.rootn R.projectn '\data\storage'];
R.data.datatype = 'NPD';
R.frqz = [4:.5:80];
R.frqzfull = [4:.5:80];
% R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'var1','var2','var2'}; % Channels in data (these will be used for comparison)
R.chsim_name = {'var1','var2','var2'}; % Channels in Sim (if hidden states exist)
R.out.tag = 'MVAR_NPD_ABC';
% Set SimAn Parameters
R.SimAn.pOptList = {'.params'};%,'.noisecov'};
R.SimAn.pOptBound = [-8 8];
R.SimAn.pOptRange = R.SimAn.pOptBound(1):.1:R.SimAn.pOptBound(2);

R.SimAn.searchN = 100;
R.SimAn.Tm = 0.8; % Initial temperature
R.SimAn.alpha = 0.98; % alpha increment
R.SimAn.rep = 64; % Repeats per temperature
R.SimAn.rtol_repeat = 0.85;
R.SimAn.rtol_converge = 0.98;
R.SimAn.ntol = 15;
R.SimAn.gradtol = [0.075 0.05];
R.SimAn.saveout = 'xobs1';
R.SimAn.maxdev = 5;
R.SimAn.jitter = 5;
R.SimAn.dSkew = 0.05;
R.SimAn.dPrec = 0.05;

% Set simulation parameters
R.IntP.intFx = @fx_simulateMVAR;
R.IntP.intFxArg ='(x,m,p,R)';
R.IntP.compFx = @compareData_100717;
R.IntP.compFxArg = '(R,sim_dat)';
R.IntP.dt = 1/200;
R.IntP.tend = 100;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);
R.IntP.Utype = 'zero'
% desired freq res:
R.obs.transFx = @constructNPDMat_190618;
R.obs.SimOrd = 8;
R.obs.brn = 0;
%(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,9,R)';

R.objfx.feattype = 'ForRev'; % 
R.objfx.specspec = 'cross'; % which part of spectra to fit
R.plot.outFeatFx = @npdplotter_110717;
R.plot.save = 'True';
R.plot.distchangeFunc = @plotDistChange_KS;
R.plot.gif.delay = 0.3;
R.plot.gif.start_t = 1;
R.plot.gif.end_t = 1;
R.plot.gif.loops = 2;




