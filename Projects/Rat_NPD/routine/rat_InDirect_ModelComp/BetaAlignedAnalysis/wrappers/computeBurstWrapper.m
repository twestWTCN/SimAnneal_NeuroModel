function [R,BB] = computeBurstWrapper(R)
load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STR->GPe','150% STR->GPe'};
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 6 8 16 18],:);
R.condcmap(6,:) = [0 0 0];

xsimMod{2} = 1e8.*xsim_HD{1}{1}; % Low strength
xsimMod{1} = 1e8.*xsim_HD{6}{1}; % Fitted
xsimMod{3} = 1e8.*xsim_HD{8}{1};% High Strength

xsimMod{5} = 1e8.*xsim_STR_GPe{1}{1}; % Low strength
xsimMod{4} = 1e8.*xsim_STR_GPe{6}{1}; % Fitted
xsimMod{6} = 1e8.*xsim_STR_GPe{10}{1};% High Strength


[R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
R.BB.thresh_prctile = 90;
BB = compute_BurstThreshold(R,BB,1,0);
R.BB.minBBlength = 1; %  Minimum burst period- cycles
BB.plot.durlogflag = 0;
BB = defineBetaEvents(R,BB);
save([R.rootn '\routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])