function [R] = BAA_sim_ConnectionSweep_v2(Rorg,modID,simtime,HD)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(Rorg,modID,simtime);

if HD == 0
    % connection list
    for CON = 1:4
        ck_1(CON,:) = [0.00001 0.125 0.25 0.5 0.75 1 1.25  1.5 3 5];
    end
    hdext = '';
elseif HD == 1
    for CON = 1:4
        ck_1(CON,:) = logspace(-1,0.7,30);
    end
    hdext = '_F1';
elseif HD == 2
    refK = logspace(-1,0.7,30);
    for CON = 1:4
        x1 = logspace(log10(refK(Rorg.betaKrange(1,CON))),log10(refK(Rorg.betaKrange(2,CON))),10);
        x2 = logspace(log10(refK(Rorg.betaKrange(2,CON))),log10(refK(Rorg.betaKrange(3,CON))),10); 

        ck_1(CON,:) = [x1(1:end-1) 1 x2(2:end)];
    end
    hdext = '_bKF';
end
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = permMod{1}.par_rep{1};
R.IntP.getNoise = 1;
[dum1,dum2,feat_sim_noise,dum3,xsim_noise] = computeSimData(R,m,uc,XBase,0);
R.IntP.getNoise = 0;

for CON = 1:4
    feat = {};
    xsim = {};
    parfor i = 1:size(ck_1,2)
        % Now Modify
        Pbase = XBase;
        if CON == 1 % Hyperdirect
            Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,i)); %
        elseif CON == 2 % Striatal-pallidal
            Pbase.A{2}(3,2) = log(exp(Pbase.A{2}(3,2))*ck_1(CON,i)); %
        elseif CON == 3 % Pallidal-subthalamo
            Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,i)); %
        elseif CON == 4 % Subthalamo-pallidal
            Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*ck_1(CON,i)); %
        end
        [r2mean,pnew,feat_sim,dum1,xsim_gl] = computeSimData(R,m,uc,Pbase,0);
        feat{i} = feat_sim;
        xsim{i} = xsim_gl;
    end
    R2 = simannealsetup_InDirect_ModelComp();
    mkdir([R2.rootn 'routine\' R2.out.oldtag '\BetaBurstAnalysis\Data'])
    
    save([R2.rootn 'routine\' R2.out.oldtag '\BetaBurstAnalysis\Data\BB_' R2.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat' hdext '.mat'],'feat')
    save([R2.rootn 'routine\' R2.out.oldtag '\BetaBurstAnalysis\Data\BB_' R2.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim' hdext '.mat'],'xsim')
    save([R2.rootn 'routine\' R2.out.oldtag '\BetaBurstAnalysis\Data\BB_' R2.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1' hdext '.mat'],'ck_1')
end

save([R2.rootn 'routine\' R2.out.oldtag '\BetaBurstAnalysis\Data\BB_' R2.out.tag '_ConnectionSweep_noise.mat'],'xsim_noise','feat_sim_noise')

