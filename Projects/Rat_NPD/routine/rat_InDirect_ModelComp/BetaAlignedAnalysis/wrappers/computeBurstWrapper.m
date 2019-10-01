function [R,BB] = computeBurstWrapper(R)
% load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim.mat'],'xsim_HD','xsim_STR_GPe')
R.BBA_path ='C:\Users\twest\Documents\Work\GitHub\BurstToolbox';
addpath( genpath(R.BBA_path));
% R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
R.condname = num2cell(1:20);
cmap = brewermap(18,'Spectral');
R.condcmap = cmap([1 4 8 16 4 18],:);
hdext = {'','_F1','_bKF','_KrSel'};

% R.condcmap(6,:) = [0 0 0];
for CON = 1:4
    for HD = 4
        if HD == 1 || HD == 2  || HD == 3
            HDM = HD;
        elseif HD == 4
            HDM = 2;
        end
        BB = [];
        R.condname = {};
        load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext{HDM} '.mat'],'ck_1')
        load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim' hdext{HDM} '.mat']); % The low density estimate
        load([R.rootn 'routine\' R.out.tag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext{HDM} '.mat']); % The low density estimate
        
        if HD == 1 || HD == 2 || HD == 3
            MList = 1:numel(xsim);
        elseif HD == 4
            MList = R.betaKrange(:,CON)';
        end
        
        xsimMod = {}; bpowr = []; fpow = []; R.condname = {};
        ip = 0;
        for i = MList
            ip = ip+1;
            xsimMod{ip} = 1e7.*xsim{i}{1};
            [bpowr(ip) b] = max(feat{i}(1,4,4,3,:));
            fpow(ip) =R.frqz(b);
            R.condname{ip} = num2str(ck_1(i),2);
        end
        
        condsel = find(bpowr<1e-8);
        
        [R,BB] = compute_BetaBursts_Simulated(R,xsimMod);
        R.BB.thresh_prctile = 75;% o85; tl 80
        BB = compute_BurstThreshold(R,BB,condsel,0);
        R.BB.minBBlength = 1; %o1 tl 1.5; %  Minimum burst period- cycles
        BB.plot.durlogflag = 0;
        if HD == 1 || HD == 3 || HD == 4
            memflag = 0;
        else
            memflag = 1; %Crop time series data (for memory)
        end
        BB = defineBetaEvents(R,BB,memflag);
        BB = getBurstStatsXConds(BB);
        BB.condlist = ck_1;
        save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BBA_' R.out.tag '_Sims_CON_' num2str(CON) hdext{HD} '.mat'],'BB','-v7.3')
    end
end