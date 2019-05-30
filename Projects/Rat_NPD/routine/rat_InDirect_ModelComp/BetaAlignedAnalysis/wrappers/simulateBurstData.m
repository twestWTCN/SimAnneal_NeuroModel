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
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN ->GPe','150% STN->GPe'};
R.condnames =  R.condname;
% R.data = feat{end};
R = BAA_sim_ConnectionSweep(R,10,1024);