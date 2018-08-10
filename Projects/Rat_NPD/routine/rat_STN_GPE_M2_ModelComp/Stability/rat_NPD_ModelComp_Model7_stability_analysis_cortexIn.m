% Model 7
% % Model 3 (M2 Feed) but with model 1 priors


%%%%%%%%%%%%%%%%%%%%%%%%
% simAnnealAddPaths()
clear ; close all
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\sim_machinery'))
% addpath(genpath('C:\Users\twest\Documents\Work\GitHub\SimAnneal_NeuroModel\Projects\Rat_NPD'))
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\TWtools\')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\bplot\')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\MEG_STN_Project')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\Neurospec\neurospec21')
% addpath('C:\spm12')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\export_fig')
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% %
% addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\sort_nat')
rng(424124)

%% Set Routine Pars
R = simannealsetup_STN_GPe_M2_ModelComp;
%% Prepare the data
R = prepareRatData_Group_NPD(R);

%% Prepare Model
[R p m uc] = MS_rat_STN_GPe_M2_ModelComp_Model7(R);

R.out.dag = 'NPD_ModComp_M7'; % 'All Cross'
R.out.tag = 'ModComp';

% load modelfit
load([R.rootn 'outputs\' R.out.tag '\NPD_ModComp_M7\modelfit_' R.out.tag '_NPD_ModComp_M7.mat'])
R = varo;
R.obs.gainmeth = R.obs.gainmeth(1);
R.obs.trans.norm = 0;
R = setSimTime(R,18);
p = R.Mfit.BPfit;
% Inverted Parameters
% uc = innovate_timeseries(R,m);
% uc{1} = uc{1}.*sqrt(R.IntP.dt);
% [xsims dum wflag J] = R.IntP.intFx(R,m.x,uc,p,m);
%
% xsims = R.obs.obsFx(xsims,m,p,R);
% [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
% R.plot.outFeatFx({R.data.feat_emp},{feat_sim},R.data.feat_xscale,R,1,[])
%
NI = 32; NJ = 32;
CswitchI = linspace(-8,8,NI); % Modulating the GPe -| STN (3,2)

Jsave= zeros(1,NI);
Jmodsave  =zeros(4,NI);
STN = zeros(1,NI);
GPe = zeros(1,NI);

j = 1;
parfor i = 1:NI
    [j i]
    pnew = p;
    pnew.A{1}(3,1) = CswitchI(i);
    
    namerz = 'none';
    rng(12312);
    % Simulate New Data
    u = innovate_timeseries(R,m);
    u{1} = u{1}.*sqrt(R.IntP.dt);
    [xsims dum wflag J] = R.IntP.intFx(R,m.x,u,pnew,m);
    eJ = eig(J{1});
    Jsave(i) =eJ(1);
    JMS = zeros(4,1);
    for mod = 1:size(m.xinds,1)+1
        if mod ~= 4
            eJ = eig(J{1}(m.xinds(mod,1):m.xinds(mod,2),m.xinds(mod,1):m.xinds(mod,2)));
            JMS(mod) = eJ(1);
        else
            eJ = eig(J{1}(m.xinds(2,1):m.xinds(3,2),m.xinds(2,1):m.xinds(3,2)));
            JMS(mod) = eJ(1);
        end
    end
    Jmodsave(:,i) = JMS;
    if wflag == 0
        try
            % Run Observer function
            xsims = R.obs.obsFx(xsims,m,pnew,R);
            % Run Data Transform
            [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
            xlim([5 8]); ylim([-5 12])
            set(gcf,'Position',[100   100   500   202])
            hAxes = gca;
            %         hAxes.XRuler.Axle.LineStyle = 'none';
            axis off
            
            feat_simjr{i} = feat_sim;
            beta_ind = find(R.data.feat_xscale > 14  & R.data.feat_xscale <24);
            STN(i) = max(squeeze(feat_sim(1,3,3,1,beta_ind)));
            GPe(i) = max(feat_sim(1,2,2,1,beta_ind));
        catch
            STN(i) = NaN;
            GPe(i) = NaN;
        end
    else
        STN(i) = NaN;
        GPe(i) = NaN;
        %             Jsave(:,i,j) = nan(1,3);
    end
    %         figure
    %         R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1,[])
end



% % close all
% % [cmap] = brewermap(128,'RdYlBu');
% % 
% % figure
% % colormap(cmap)
% % subplot(2,1,1)
% % colGridCon(CswitchI,CswitchJ,(GPe),4)
% % xlabel('M2 \rightarrow  STN'); ylabel('GPe \rightarrow STN')
% % title('GPe Beta Power','FontSize',14);% xlim([-2 4]);ylim([-3 1]); caxis([-2 1.5])
% % a = gca;
% % a.FontSize = 16;
% % 
% % subplot(2,1,2)
% % colGridCon(CswitchI,CswitchJ,(STN),4)
% % xlabel('M2 \rightarrow  STN'); ylabel('GPe \rightarrow STN')
% % title('STN Beta Power','FontSize',14);% xlim([-2 4]);ylim([-3 1]); caxis([-2 1.5])
% % a = gca;
% % a.FontSize = 16;
% % set(gcf,'Position',[1286         129         514         858])
% % 
% % 
% % figure
% % colormap(cmap)
% % subplot(2,1,1)
% % colGridCon(CswitchI,CswitchJ,(squeeze(real(Jsave(1,:,:)))),-1.5:0.5:0.5)
% % xlabel('M2 \rightarrow  STN'); ylabel('GPe \rightarrow STN')
% % title('Real Part','FontSize',14);% xlim([-2 4]);ylim([-3 1]);
% % caxis([-2 1])
% % a = gca;
% % a.FontSize = 16;
% % set(gcf,'Position',[1286         129         514         858])
% % 
% % subplot(2,1,2)
% % colGridCon(CswitchI,CswitchJ,(squeeze(imag(Jsave(1,:,:)))),0:0.25:1)
% % xlabel('M2 \rightarrow  STN'); ylabel('GPe \rightarrow STN')
% % title('Imaginary EV','FontSize',14);% xlim([-2 4]);ylim([-3 1]); caxis([-2 1.5])
% % a = gca;
% % a.FontSize = 16;
% % set(gcf,'Position',[1286         129         514         858])
