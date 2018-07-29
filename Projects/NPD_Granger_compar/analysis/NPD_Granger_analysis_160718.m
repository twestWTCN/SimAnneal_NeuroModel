clear ; close all
R = simannealsetup_comparNPD_Granger_170118;
R.out.tag = 'comparNPD_Granger';
R.out.dag = '100718';
R.out.dag = [R.out.dag 'cross'];
% R.SimAn.pOptList = {'.int{src}.T','.int{src}.G','.int{src}.S','.S','.C','.obs.LF','.A','.D'}; %,'.S','.obs.LF'}; ,'.int{src}.S','.S' %,'.C','.obs.LF'}; % ,'.obs.mixing','.C','.D',

load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelspec_' R.out.tag '_' R.out.dag '.mat'])
m = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\modelfit_' R.out.tag '_' R.out.dag '.mat'])
R = varo;
load([R.rootn 'outputs\' R.out.tag '\' R.out.dag '\parBank_' R.out.tag '_' R.out.dag '.mat'])
parBank = varo;

% Need m,p, parBank
p = R.Mfit.Pfit;

addpath(genpath('C:\Users\twest\Documents\Work\Github\SimAnneal_NeuroModel\sim_machinery'))
% Create parameter index map
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form


% Re-compute optimised parbank
R.analysis.modEvi.N = 100;
R.analysis.modEvi.eps = 0.70;
parOptBank = parBank(:,parBank(end,:)>=R.analysis.modEvi.eps);
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
triangle_plot_mdist_3node(R,p,m,xf)

%% Compute Model Evidence

figure(2)
modelProbs(m.x,m,p,R,R.out.dag)

%% Plot Simulated Features with Confidence intervals
figure(3)
PlotFeatureConfInt_gen(R,R.out.dag)
saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '1\featConfInts'],'-jpg',1,'-r200',2);
optP = getOptParMean(m,p,R,R.out.dag);



%% Plot Model Output
Plot_Connectivity_Effects


% 
% 
% %% Parameter Sweep across STN/GPe Subcircuit Connections
% parsweep.R = '.A{1}(2,1)';
% parsweep.Rname = '1->2';
% parsweep.Q = '.A{1}(3,1)';
% parsweep.Qname = '1->3';
% 
% % Conduct sweep
% parsweep = modelBetaParSweep(m,optP,parsweep,R);
% pathstr = [R.rootn 'analysis\parsweeps\'];
% mkdir(pathstr);
% save([pathstr '\A_STNGPe_AGPeSTN_parsweep'],'parsweep');
% 
% % Plot heatmap
% figure(4)
% parSweepPlot(R,parsweep)
%  saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '\parSweepSTNGPe.jpg'],'-jpg',1,'-r200',2);
