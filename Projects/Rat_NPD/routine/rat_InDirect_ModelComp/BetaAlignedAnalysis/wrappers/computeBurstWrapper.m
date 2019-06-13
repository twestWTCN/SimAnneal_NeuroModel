function [R,BB] = computeBurstWrapper(R)
load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
% R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
R.condname = num2cell(1:20)
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 4 8 16 4 18],:);
% R.condcmap(6,:) = [0 0 0];

% xsimMod{2} = 1e7.*xsim_HD{1}{1}; % Low strength (uV)
% xsimMod{1} = 1e7.*xsim_HD{6}{1}; % Fitted
% xsimMod{3} = 1e7.*xsim_HD{10}{1};% High Strength
% 
% xsimMod{5} = 1e7.*xsim_STR_GPe{5}{1}; % Low strength
% xsimMod{4} = 1e7.*xsim_STR_GPe{6}{1}; % Fitted
% xsimMod{6} = 1e7.*xsim_STR_GPe{10}{1};% High Strength

for i = 1:10
xsimMod{i} = 1e7.*xsim_HD{i}{1};
end
for i = 11:20
xsimMod{i} = 1e7.*xsim_STR_GPe{i-10}{1};
end

[R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
R.BB.thresh_prctile = 90;% o85; tl 80
BB = compute_BurstThreshold(R,BB,[1:20],0);
R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
BB.plot.durlogflag = 0;
BB = defineBetaEvents(R,BB);
save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims.mat'])
