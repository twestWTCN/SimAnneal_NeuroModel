clear ; close all
R = simannealsetup_NPD_300817_OFF;

d = 'trialanalysisOFF';
% Need m,p, parBank
% load('emerg_save.mat','p')

% load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\M_OFF_FRESH.mat')
% 
% R.out.tag = 'NPD_Final_JNPPaper_lesion_fresh';
% R.out.dag = '2018627';
% 
% load([R.rootn '\' R.projectn '\outputs\' R.out.tag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
% R = varo;
% R.out.dag = '2018627';
% 
% load([R.rootn '\' R.projectn '\outputs\' R.out.tag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
% parBank = varo;
% clear varo 

load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\emerg_save.mat')
p_off = xobs1.Mfit.Pfit;
R.obs.logdetrend = 1; % Detrend Autospectra



fx = R.data.feat_emp;
fxa(1,:,:,:,:) = fx;
R.data.feat_emp = fxa;
R.condnames = {'off'}
m.outstates = {[0 0 0 0 0 0 1 0]  [1 0]  [1 0]  [1 0]  [1 0]  [1 0]};
R.obs.outstates = find([m.outstates{:}]);
for i=1:numel(R.chloc_name)
    R.obs.obsstates(i) = find(strcmp(R.chloc_name{i},R.chsim_name));
end
addpath(genpath('C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery'))
pathstr = [R.rootn '\outputs\' R.out.tag '\' d '\OptSaves\'];
% load optimised model
% load([pathstr 'modelfit_' R.out.tag '_' d '.mat'])

% Create parameter index map
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form


% Re-compute optimised parbank
R.analysis.modEvi.N = 100;
R.analysis.modEvi.eps = 0.09;
parOptBank = parBank(:,parBank(end,:)>-0.1);
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

figure(2)
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
