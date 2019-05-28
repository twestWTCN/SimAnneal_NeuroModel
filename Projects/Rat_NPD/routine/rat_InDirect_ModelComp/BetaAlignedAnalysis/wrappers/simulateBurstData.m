function R = simulateBurstData(R)
%% Simulate Data
R.bandinits = {'\alpha','\beta_1','\beta_2'};
R.condcmap = linspecer(3);
R.cohband = 2;
R.obs.trans.norm = 0;
% R = BAA_sim_discreteModels(R,10,250);
% OR load
% load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],'xsimMod','feat')

% Setup condition variables
R.condname = {'Full','1% M2->STN','150% M2->STN','1% STR->GPe','150% STR->GPe','Empirical'};
R.condnames =  R.condname;
% R.data = feat{end};
R = BAA_sim_ConnectionSweep(R,10,500);