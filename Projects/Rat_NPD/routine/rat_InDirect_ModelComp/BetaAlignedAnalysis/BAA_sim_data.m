R.out.tag = 'InDrt_ModComp';
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
%% get data from model
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,permMod,xsimMod] = getSimData_v2(R,8,25);

clear r2mean feat_sim xsims
Alist = -4:0.5:3;
for i = 1:length(Alist)
    pnew = permMod{1}.par_rep{1};
    pnew.A{1}(4,1) = Alist(i);
    [r2mean{i},dum,feat_sim{i},xsims{i}] = computeSimData(R,m,pnew,30,1); 
    fx = squeeze(feat_sim{i}(1,4,4,1,:));
    bp(i) = sum(fx(R.frqz>14 & R.frqz<24));
end  
    
R.out.tagOld = 'rat_InDirect_ModelComp';
mkdir([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data'])
save([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],...
    'permMod','xsimMod')
% OR
load([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'])
%% Get autospectra (for frequency finding)
AS{1} = squeeze(permMod{1}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
AS{2} = squeeze(permMod{2}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
% AS{3} = squeeze(permMod{3}.feat_rep{1}(1,:,4,1,:)); % get STN autospectra
clear('permMod','permModHD','permModSTR')

