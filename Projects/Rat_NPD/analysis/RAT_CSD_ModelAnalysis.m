clear ; close all
R = simannealsetup_170817;
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\crossfit_delays_pmuparBankparBankOpts.mat')
load([R.rootn 'data\storage\datafeat_npd']);
R.data.feat_emp = meannpd_data;
R.data.feat_xscale = F_data;


%% Compute Model Evidence
R.Mfit = R_out.Mfit;
R.out.tag = 'CSD_ABC_neatmodel_iterate';

R.analysis.modEvi.N = 1000;
R.analysis.modEvi.eps = -0.1;
d =  sprintf('%d',[R_out.d(1:3)]);
R.parOptBank = parOptBank;
figure(1)
modelProbs(m.x,m,p,R,d)

% %% Plot Fit with Confidence intervals
figure(2)
PlotFeatureConfInt(R,d)
saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '2\featConfInts'],'-jpg',1,'-r200',2);

optP = getOptParMean(m,p,R,d);

parsweep.R = '.A{1}(3,4)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Q = '.A{2}(4,3)';
parsweep.Qname = 'GPe-|STN Connections Strength';

parsweep = modelBetaParSweep(m,optP,parsweep,R);
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STNGPe_AGPeSTN_parsweep'],'parsweep');
figure(2)
parSweepPlot(R,parsweep)
saveallfiguresFIL_n([R.rootn '\analysis\' R.out.tag '\parSweepSTNGPe.jpg'],'-jpg',1,'-r200',2);
