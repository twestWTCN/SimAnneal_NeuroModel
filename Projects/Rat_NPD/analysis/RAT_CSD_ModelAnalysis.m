clear ; close all
R = simannealsetup_170817;
addpath(genpath('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\sim_machinery'))
load('C:\Users\Tim\Documents\Work\GIT\SimAnneal_NeuroModel\Projects\Rat_NPD\Saves\cross_fit_pars_Rout_p_m_u.mat')
load([R.rootn 'data\storage\datafeat_npd']);
R.data.feat_emp = meannpd_data;
R.data.feat_xscale = F_data;


%% Compute Model Evidence
R.Mfit = R_out.Mfit
R.out.tag = 'CSD_ABC_neatmodel1';

R.analysis.modEvi.N = 250;
R.analysis.modEvi.eps = -0.5;
d = '2017818';
% modelProbs(m.x,m,p,R,d)

% %% Plot Fit with Confidence intervals
% PlotFeatureConfInt(R,d)


optP = getOptParMean(m,p,R,d);

parsweep.R = '.A{1}(3,4)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Q = '.A{2}(4,3)';
parsweep.Qname = 'GPe-|STN Connections Strength';

parsweep = modelBetaParSweep(m,optP,parsweep,R);
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STNGPe_AGPeSTN_parsweep'],'parsweep');


