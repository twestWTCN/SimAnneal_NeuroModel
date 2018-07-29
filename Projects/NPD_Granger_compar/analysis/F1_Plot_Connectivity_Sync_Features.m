%% Plot Model Output
clear
load('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\NPD_Granger_compar\analysis\analysis_save_figure1_connectioneffects.mat')
close all
R = setSimTime(R,200);
R.obs.gainmeth = R.obs.gainmeth(1:2);
% Set up for 2 nodes
R.obs.LF = R.obs.LF(1:2);
R.obs.obsstates = [1 2];
R.chloc_name = R.chloc_name(1:2);
R.chsim_name = R.chsim_name(1:2);

m.m = 2;
m.x = m.x(1:2);
m.Gint = m.Gint(1:2);
m.Tint = m.Tint(1:2);
m.dipfit.model = m.dipfit.model(1:2);
m.outstates = m.outstates(1:2);
R.obs.outstates = R.obs.outstates(1:2)
m.xinds = m.xinds(1:2,:);
m.uset.p.covar = m.uset.p.covar(1:2,1:2);
m.n = sum(numel([m.outstates{:}]));
p.A{1} = p.A{1}(1:2,1:2);p.A{2} = p.A{2}(1:2,1:2);
p.C = [0 0]';
p.obs.LF = [0 0];
p.D = zeros(2);
p.int = p.int(1:2);
p.int{2} = p.int{1};
% Nought Connectivity
Anought = repmat(-32,2);
Acon = Anought;
Cswitch = 2.6; %linspace(0,6,64);
% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);
for jr = 1:2
    for i = 1:length(Cswitch)
        pnew= p;
        Acon = Anought;
        if jr == 2
            Acon(2,1) = Cswitch(i);
            %     elseif jr==2
            %         Acon(1,2) = Cswitch(i);
        elseif jr==1
            Acon(2,1) = Cswitch(i);
            Acon(1,2) = Cswitch(i);
        end
        
        pnew.A{1} = Acon;
        pnew.A{2} = Anought;
        namerz = 'none';
        rng(124);
        % Simulate New Data
        u = innovate_timeseries(R,m);
        u{1} = u{1}.*sqrt(R.IntP.dt);
        xsims = R.IntP.intFx(R,m.x,u,pnew,m);
        % Run Observer function
        if isfield(R.obs,'obsFx')
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        
        figure(1)
        for ip =1:2
            plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}(ip,:)'+((jr-1)*5),'color',cmap(ip,:)); hold on
        end
        xlim([25 28]); xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
        % Run Data Transform
        if isfield(R.obs,'transFx')
            [~,feat_sim{jr}] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,11,R);
        else
            feat_sim = xsims; % else take raw time series
        end
    end
end
figure(3)
R.plot.outFeatFx({feat_sim{1}},{feat_sim{2}},R.data.feat_xscale,R,1,[])
subplot(2,2,2);legend({'X1 <-> X2 (A.S.)','X1 <-> X2','X1 -> X2 (A.S.)','X1 -> X2','X2 <- X1 (A.S.)','X2 <- X1'});


figure(1)
hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none';
axis off
set(gcf,'Position',[337         673        1449         331])
% figure(2)
% hAxes = gca;
% hAxes.XRuler.Axle.LineStyle = 'none';
% axis off
% set(gcf,'Position',[337         673        1449         331])

% legend({'X1