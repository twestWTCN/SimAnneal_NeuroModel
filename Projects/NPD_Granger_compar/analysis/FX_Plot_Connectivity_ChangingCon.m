%% Plot Model Output
clear
load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\NPD_Granger_compar\analysis\analysis_save_figure1_connectioneffects.mat')
close all
R = setSimTime(R,32);
for i = 1:4
pnew =optP;
% pnew.int{2} =pnew.int{1};
% pnew.int{3} =pnew.int{1};
if i ==1
    namerz = 'fitted';
elseif i == 2
    At = pnew.A{1};
    At(find(tril(At,-1))) = -32;
    pnew.A{1} = At;
    namerz = 'no_forward';
elseif i==3
    At = pnew.A{1};
    At(find(triu(At,1))) = -32;
    pnew.A{1} = At;
    namerz = 'no_reverse';
elseif i==4
    pnew.A{1} = pnew.A{1}-3; %repmat(-32,3,3);
    pnew.int{2} =pnew.int{1};
    pnew.int{3} =pnew.int{1};
    namerz = 'none';
end


% Simulate New Data
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
xsims = R.IntP.intFx(R,m.x,u,pnew,m);
% Run Observer function
if isfield(R.obs,'obsFx')
    xsims = R.obs.obsFx(xsims,m,pnew,R);
end
xsimsave{i} = xsims;
% Run Data Transform
if isfield(R.obs,'transFx')
    [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
else
    feat_sim = xsims; % else take raw time series
end
% Compare Pseudodata with Real
r2mean  = R.IntP.compFx(R,feat_sim);
figure
plot(xsims{1}')
figure
 R.plot.outFeatFx({R.data.feat_emp},{feat_sim},R.data.feat_xscale,R,1,[])
% G = gcf;
% G.Name = namerz;
end
