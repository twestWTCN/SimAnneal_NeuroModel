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
Cswitch = linspace(0,6,64);
% plotting ops
Lss = {'-','--'}; %,'-.'};
cmap = linspecer(3);
for jr = 1:2
    parfor i = 1:length(Cswitch)
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
        
        figure
        plot(linspace(0,length(xsims{1})*R.IntP.dt,length(xsims{1})),xsims{1}')
        % Run Data Transform
        if isfield(R.obs,'transFx')
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
        else
            feat_sim = xsims; % else take raw time series
        end
        %          R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
        
        % NPD
        for  nd = 1:3
            npdmax(nd,i) = max(feat_sim(1,1,2,nd,:));
            %         npdmax(2,nd,i) = max(feat_sim(1,2,1,nd,:));
        end
        
        % DFA
        DFAP= [];
        DFAP(1) = 1/R.IntP.dt; DFAP(2) = (8/30); DFAP(3) = (length(xsims{1})/8)*R.IntP.dt; DFAP(4) = 50; DFAP(5) =0;
        for ns = 1:2
            [bmod win evi alpha] = peb_dfa_TW_2018(abs(hilbert(xsims{1}(ns,:))),DFAP,-4,0);
            DFAevi(ns,i) = evi;
            DFAalpha(ns,i) = alpha;
            stdAbs(ns,i) = std(abs(hilbert(xsims{1}(ns,:))));
        end
        
        % DFAPS
        DFAP= [];
        DFAP(1) = 1/R.IntP.dt; DFAP(2) = (8/30); DFAP(3) = (length(xsims{1})/8)*R.IntP.dt; DFAP(4) = 50; DFAP(5) =0;
        x1 = angle(hilbert(xsims{1}(1,:))); x2 = angle(hilbert(xsims{1}(2,:)));
        x = diff(unwrap(x1-x2));
        [bmod win evi alpha] = peb_dfa_TW_2018(x,DFAP,-4,0);
        DFAPSevi(i) = evi;
        DFAPSalpha(i) = alpha;
        stdDPS(i) = std(x);
        disp(i)
    end
    
    %% NPD
    figure(1)
    X = squeeze(npdmax);
    for nd = 1:3
        plot(Cswitch,X(nd,:),'color',cmap(nd,:),'LineWidth',2,'LineStyle',Lss{jr});
        hold on
    end
    %% DFA
    figure(2)
    subplot(2,2,1)
    for nd = 1:2
        X = squeeze(DFAevi);
        plot(Cswitch,X(nd,:),'color',cmap(nd,:),'LineWidth',2,'LineStyle',Lss{jr});
        hold on
    end
    subplot(2,2,3)
    for nd = 1:2
        X = squeeze(DFAalpha);
        plot(Cswitch,X(nd,:),'color',cmap(nd,:),'LineWidth',2,'LineStyle',Lss{jr});
        hold on
    end
    
    %% DFAPS
    figure(2)
    subplot(2,2,2)
    X = squeeze(DFAPSevi);
    plot(Cswitch,X(1,:),'color',cmap(nd,:),'LineWidth',2,'LineStyle',Lss{jr});
    hold on
    subplot(2,2,4)
    X = squeeze(DFAPSalpha);
    plot(Cswitch,X(1,:),'color',cmap(nd,:),'LineWidth',2,'LineStyle',Lss{jr});
    hold on
    
     %% NPD
    figure(3)
    X = squeeze(npdmax);
    for nd = 1:3
        plot(Cswitch,stdDPS,'color',cmap(3,:),'LineWidth',2,'LineStyle',Lss{jr});
        hold on
        plot(Cswitch,stdAbs(1,:),'color',cmap(1,:),'LineWidth',2,'LineStyle',Lss{jr});
        plot(Cswitch,stdAbs(2,:),'color',cmap(2,:),'LineWidth',2,'LineStyle',Lss{jr});
    end   
end
set(gca, 'YScale', 'log')
figure(1)
legend({'X1 <-> X2','X1 -> X2','X2 -> X1','X1 <-> X2 (A.S.)','X1 -> X2 (A.S.)','X2 -> X1 (A.S.)'});
xlabel('Coupling'); ylabel('NPD')
set(gcf,'Position',[288   408   701   618])

figure(2)
subplot(2,2,1)
xlabel('Coupling'); ylabel('DFA Evidence');xlim([0 6])
subplot(2,2,3)
xlabel('Coupling'); ylabel('DFA Alpha');xlim([0 6])
subplot(2,2,2)
xlabel('Coupling'); ylabel('DFA Evidence');xlim([0 6])
subplot(2,2,4)
xlabel('Coupling'); ylabel('DFA-PS Alpha');xlim([0 6])
set(gcf,'Position',[288   408   701   618])

figure(3)
legend({'std dRP(X1,X2)','std |H(X1)|','std |H(X2)|'});
xlabel('Coupling'); ylabel('log std')
set(gcf,'Position',[288   408   701   618])

