function [R] = BAA_sim_ConnectionSweep(R,modID,simtime,HD)
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);

if HD == 0
    % connection list
    ck_1 = [0.00001 0.1 0.25 0.5 0.75 1 1.5  2  2.1 2.5];
%     ck_2 = [0.00001 0.1 0.25 0.5 0.75 1 1.05 1.10 1.20 1.3]; % STN - > GPe
    ck_2 = [0.00001 0.1 0.25 0.5 0.75 1 1.025  1.1  1.125 1.5]; % GPe - > STN
    hdext = '';
elseif HD == 1
    ck_1 = logspace(-1,0.6,30);
%     ck_2 = logspace(-1,0.14,30); % STN - > GPe
    ck_2 = logspace(-0.5,0.15,30); % GPe - > STN
    hdext = '_F1';
end
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = permMod{1}.par_rep{1};
parfor i = 1:numel(ck_1)
    % Now Modify
    % Model 1 (hyperdirect)
    Pbase = XBase;
    Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(i)); %
    [r2mean,pnew,feat_sim,dum,xsim_gl] = computeSimData(R,m,uc,Pbase,0);
    feat_HD{i} = feat_sim;
    xsim_HD{i} = xsim_gl;
    % Model 2 (striatopallidal)
    Pbase = XBase;
%     Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*ck_2(i)); %
Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_2(i)); %
    
    [r2mean,pnew,feat_sim,dum,xsim_gl] = computeSimData(R,m,uc,Pbase,0);
    feat_STR_GPe{i} = feat_sim;
    xsim_STR_GPe{i} = xsim_gl;
end

R = simannealsetup_InDirect_ModelComp();
mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_feat' hdext '.mat'],'feat_HD','feat_STR_GPe')
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep_xsim' hdext '.mat'],'xsim_HD','xsim_STR_GPe')