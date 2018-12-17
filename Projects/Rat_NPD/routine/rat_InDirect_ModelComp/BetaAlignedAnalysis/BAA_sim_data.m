function R = BAA_sim_data(R,modID,simtime)
[R,permMod,xsimMod] = getSimModelData_v2(R,modID,simtime); 

% Now you have posterior expectations you can use the inside computeSim
% Data Function
% To run through changing the connection strengths:
% clear r2mean feat_sim xsims
% Alist = -4:0.5:3;
% for i = 1:length(Alist)
%     pnew = permMod{1}.par_rep{1};
%     pnew.A{1}(4,1) = Alist(i); %hyperdirect input
%     [r2mean{i},dum,feat_sim{i},xsims{i}] = computeSimData(R,m,pnew,30,1); 
%     fx = squeeze(feat_sim{i}(1,4,4,1,:));
%     bp(i) = sum(fx(R.frqz>14 & R.frqz<24));
% end  
    
R.out.tagOld = 'rat_InDirect_ModelComp';
mkdir([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data'])
save([R.rootn '\routine\' R.out.tagOld '\BetaBurstAnalysis\Data\BB_' R.out.tag '_Sims.mat'],...
    'permMod','xsimMod')
