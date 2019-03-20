function R = PB_scheme_setup()
R.projectn = 'ABC_PBSchema';
switch getenv('computername')
    case 'SFLAP-2'
        usname = 'Tim'; gitpath = '\Documents\Work\GIT'; madpath = 'MATLAB_ADDONS';
        spmpath = 'C:\Users\Tim\Documents\spm12';
    case 'FREE'
        usname = 'twest'; gitpath = '\Documents\Work\GitHub'; madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\spm12';
    case 'DESKTOP-94CEG1L'
        usname = 'timot'; gitpath =  '\Documents\GitHub'; madpath = 'Work\MATLAB ADDONS';
        spmpath = 'C:\Users\timot\Documents\GitHub\spm12';
end
R.rootn = ['C:\Users\' usname gitpath '\SimAnneal_NeuroModel\Projects\' R.projectn '\'];

R.out.tag = 'ABC_PBSchema';
R.out.dag = 'PB_scheme1';
%% Setup Simulation parameters
R.condnames = 'STIM';

% Optimization
R.SimAn.searchMax = 200;
R.SimAn.rep = 132;
R.SimAn.pOptList = {'.IWS','.SCRate','.EPSP_amp','.EPSP_ampJit'}; %'.EPSP_Tdecay'
R.SimAn.jitter = 2;
R.SimAn.pOptBound = [-8 8];
R.SimAn.pOptRange = -3:.01:3;
R.SimAn.convIt = 0.0001;
% Simulation
R.IntP.intFx = @PB_schema_simulate;
R.IntP.compFx = @PB_compareData;
R.IntP.Utype = 'zero';
 
R.IntP.dt = 0.0001; % Not needed!!
R.frqzfull = 500; % Not needed!!
R.IntP.nt = 1e3; % Not needed!!
% Setup Setup
R.model.nPulses = 5000;

% Time Scale
R.model.fsamp = 4048;
R.model.tend = 500;
R.model.t_in = linspace(0,R.model.tend,(R.model.fsamp*R.model.tend));

% Data tranform
R.obs.transFx = @PB_DataTransform;

R.obs.glist = 1;

R.trans.pirange = linspace(-pi,pi,10);
R.trans.amprange = 0.1:0.2:3.2;
% Plot Outputs
R.plot.outFeatFx = @PB_PlotOutputs;

