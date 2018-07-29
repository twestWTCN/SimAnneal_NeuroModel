clear ; close all
R = simannealsetup_NPD_300817_OFF;
addpath(genpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\boundedline-pkg'))
d = 'trialanalysisOFF';
% Need m,p, parBank
% load('emerg_save.mat','p')

% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\M_OFF_FRESH.mat')
% 
R.out.tag = 'NPD_Final_JNPPaper';
R.out.dag = '180718_COND1_auto';

load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;
p = R.Mfit.Pfit;
R.obs.logdetrend = 1; % Detrend Autospectra


% Create parameter index map
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form


% Re-compute optimised parbank
R.analysis.modEvi.N = 100;
R.analysis.modEvi.eps = 0.09;
parOptBank = parBank(:,parBank(end,:)>-0.4);
R.parOptBank = parOptBank;

for i = 1:size(pIndMap,1)
    x = R.parOptBank(pIndMap(i),:); % choose row of parameter values
    copU(i,:) = ksdensity(x,x,'function','cdf'); % KS density estimate per parameter
    xf(i,:) = x;
end

[Rho,nu] = copulafit('t',copU','Method','ApproximateML'); % Fit copula
% Save outputs that specify the copula
R.Mfit.Rho = Rho;
R.Mfit.xf = xf;
R.Mfit.nu = nu;


%% Triangle plot of MV posterior distribution
figure(1)
triangle_plot_mdist(R,p,m,xf)

%% Compute Model Evidence

figure(2);
R.analysis.modEvi.N = 500;
modelProbs(m.x,m,p,R,d)

%% Plot Simulated Features with Confidence intervals
figure(3)
PlotFeatureConfInt_gen(R,d)
saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '1\featConfInts'],'-jpg',1,'-r200',2);

optP = getOptParMean(m,p,R,d);

%% Parameter Sweep across STN/GPe Subcircuit Connections
parsweep.R = '.A{1}(3,4)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Q = '.A{2}(4,3)';
parsweep.Qname = 'GPe-|STN Connections Strength';

% Conduct sweep
parsweep = modelBetaParSweep(m,optP,parsweep,R);
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STNGPe_AGPeSTN_parsweep'],'parsweep');

% Plot heatmap
figure(4)
parSweepPlot(R,parsweep)
 saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '\parSweepSTNGPe.jpg'],'-jpg',1,'-r200',2);
