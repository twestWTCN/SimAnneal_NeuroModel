function [R] = BAA_sim_ConnectionSweep(R,modID,simtime)
modID = 12; simtime = 32;
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);

% connection list
ck = [0.01 0.1 0.25 0.5 0.75 1 1.5 2 4 10 100];

parfor i = 1:numel(ck);
    % Now Modify
    % Model 1 (hyperdirect)
    Pbase = permMod{1}.par_rep{1};
    Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck(i)); % 
    [r2mean,pnew,feat_sim] = computeSimData(R,m,Pbase,0);
    feat_HD{i} = feat_sim;
    
    % Model 2 (striatopallidal)
    Pbase = permMod{1}.par_rep{1};
    Pbase.A{1}(4,1) = log(exp(Pbase.A{2}(3,2))*ck(i)); % 
    [r2mean,pnew,feat_sim] = computeSimData(R,m,Pbase,0);
    feat_STR_GPe{i} = feat_sim;    
end

R = simannealsetup_InDirect_ModelComp();
mkdir([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data'])
save([R.rootn 'routine\' R.out.oldtag '\BetaBurstAnalysis\Data\BB_' R.out.tag '_ConnectionSweep.mat'],'feat_HD','feat_STR_GPe')
