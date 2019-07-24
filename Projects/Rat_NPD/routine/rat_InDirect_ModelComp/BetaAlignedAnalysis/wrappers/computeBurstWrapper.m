function [R,BB] = computeBurstWrapper(R)
% load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
% R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
R.condname = num2cell(1:20);
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 4 8 16 4 18],:);
% R.condcmap(6,:) = [0 0 0];
for CON = 1:4
    ck_1 = [0.00001 0.125 0.25 0.5 0.75 1 1.5  2 4 5];
    BB = []; xsimMod = {};
    R.condname = {};
    
    load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim.mat']); % The low density estimate

    for i = 1:10
        xsimMod{i} = 1e7.*xsim{i}{1};
            R.condname{i} = num2str(ck_1(i),2);
    end
    
    [R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
    R.BB.thresh_prctile = 90;% o85; tl 80
    BB = compute_BurstThreshold(R,BB,[1:10],0);
    R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
    BB.plot.durlogflag = 0;
    BB = defineBetaEvents(R,BB);
    
    save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) '.mat'])
end