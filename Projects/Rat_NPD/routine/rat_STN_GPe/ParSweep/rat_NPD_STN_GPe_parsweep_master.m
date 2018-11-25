clear ; close all

%  simAnnealAddPaths()
rng(6453)

%% Set Routine Pars
R = simannealsetup_InDirect_ModelComp;

modID = 1;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
R.analysis.modEvi.eps = prctile(parBank(end,:),75); % Get models above 75th percentile
R.analysis.modEvi.N = 1000;
optP = getOptParMean(m,p,R,parBank);

%% Get Sample Data
[R,permMod,xsimMod] = getSimData_v2(R,1:3,25);
cmap = linspecer(2);
figure(1)
set(gcf,'Position',[680   301   337   797])
for i = 1:3
    subplot(3,1,i)
    p = plot(R.IntP.tvec_obs,xsimMod{i}{1}{1}'+[0 6],'LineWidth',1.5);
    hold on
    plot([10.1 10.6],[10 10],'k-','LineWidth',2)
    p(1).Color = cmap(1,:);
    p(2).Color = cmap(2,:);
    xlim([10 12]); ylim([-5 12])
end
%% Set expected values
% GPE Prior means
GPe.G = [2]*200;   % synaptic connection strengths
GPe.T = [8];       % synaptic time constants;
GPe.E  = [.2 .2 -.2 -.2]*10000;             % GPE connections

% STN Prior means
STN.G = [2]*200;   % synaptic connection strengths
STN.T = [4];       % synaptic time constants ;
STN.E  = [.2 .2 -.2 -.2]*10000;             % STN connections
contmap = [5 8.5 10 15 20 30; % hz contours
    -3  -2 -1 -0.75 -0.5 -0.1];  % log(amp) contours

%% Parameter Sweep across STN/GPe Subcircuit Connections
N = 32;
parsweep.N = N;
parsweep.R = '.A{1}(1,2)';
parsweep.Rname = 'STN-> GPe Connections Strength (kHz)';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Rlist_tick = STN.E(1)*exp(parsweep.Rlist);
parsweep.Q = '.A{2}(2,1)';
parsweep.Qname = 'GPe-|STN Connections Strength (kHz)';
parsweep.Qlist = linspace(-2,3,N);
parsweep.Qlist_tick = GPe.E(1)*exp(parsweep.Qlist);
parsweep.InvertXY = [optP.A{1}(1,2) optP.A{2}(2,1)];

% Conduct sweep
R.obs.gainmeth = {'unitvar'};
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R = simannealsetup_NPD_STN_GPe;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
pathstr = [R.rootn 'outputs\' R.out.tag '\' R.out.dag];
mkdir(pathstr);
save([pathstr '\'  R.out.dag '_A_parsweep'],'parsweep');




%% Parameter Sweep Time CoOnstants
% Reload Data
R = simannealsetup_NPD_STN_GPe;
modID = 1;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
R.analysis.modEvi.eps = -0.5;
R.analysis.modEvi.N = 1000;
optP = getOptParMean(m,p,R,parBank);

N = 32;
parsweep.N = N;
parsweep.R = '.int{2}.T';
parsweep.Rname = 'STN Time Constant (ms)';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Rlist_tick = STN.T.*exp(parsweep.Rlist);
parsweep.Q = '.int{1}.T';
parsweep.Qname = 'GPe Time Constant (ms)';
parsweep.Qlist = linspace(-2,3,N);
parsweep.Qlist_tick = GPe.T.*exp(parsweep.Rlist);
parsweep.InvertXY = [optP.int{2}.T optP.int{1}.T];
% Conduct sweep
R.obs.gainmeth = {'unitvar'};
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
pathstr = [R.rootn 'outputs\' R.out.tag '\' R.out.dag];
mkdir(pathstr);
save([pathstr '\'  R.out.dag '_STN_GPe_parsweep'],'parsweep');

%% NOW Plot heatmap
load([pathstr '\'  R.out.dag '_A_parsweep'],'parsweep');
figure(4)
set(gcf,'Position',[680   319   1000   1000])
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
        parsweep.contmap = contmap(1,:);
    elseif i == 2
        parsweep.plot.feat = log10(squeeze(parsweep.frqPowBank(2,:,:)));
        parsweep.contmap = contmap(2,:);
    end
    subplot(2,2,i)
    parSweepPlot(R,parsweep)
    
    g = gca;
    g.XTickLabel = sprintfc('%.1f',(STN.E(1).*exp(g.XTick))./1000);
    g.YTickLabel = sprintfc('%.1f',(GPe.E(1).*exp(g.YTick))./1000);
    
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    if i==1
        caxis([5 28])
        c.Label.String = 'Frequency (Hz)';
        title({['Effect of Connection Strength Upon ']; [R.chsim_name{2} ' Frequency']},'FontSize',10,'fontweight','bold')
    elseif i==2
        caxis([-3 0.1])
        c.Label.String = 'log_{10}( Band Power )';
        title({['Effect of Connection Strength Upon '];[R.chsim_name{2} ' Power']},'FontSize',10,'fontweight','bold')
    end
end

load([pathstr '\'  R.out.dag '_STN_GPe_parsweep'],'parsweep');
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
        parsweep.contmap = contmap(1,:);
    elseif i == 2
        parsweep.plot.feat = log10(squeeze(parsweep.frqPowBank(2,:,:)));
        parsweep.contmap = contmap(2,:);
    end
    subplot(2,2,i+2)
    parSweepPlot(R,parsweep)
    
    g = gca;
    g.XTickLabel = sprintfc('%.1f',(STN.T(1).*exp(g.XTick)));
    g.YTickLabel = sprintfc('%.1f',(GPe.T(1).*exp(g.YTick)));
    
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    if i==1
        caxis([5 28])
        c.Label.String = 'Frequency (Hz)';
        title({['Effect of Time Constants Upon '];[R.chsim_name{2} ' Frequency']},'FontSize',10,'fontweight','bold')
    elseif i==2
        caxis([-3 0.1])
        c.Label.String = 'log_{10}( Band Power )';
        title({['Effect of Time Constants Upon '];[R.chsim_name{2} ' Power']},'FontSize',10,'fontweight','bold')
    end
end
% set(gcf,'Position',[610.0000  157.5000  950.0000  635.5000])
% annotation(gcf,'textbox',...
%     [0.152052631578947 0.962399596025029 0.745842105263159 0.0299837925445702],...
%     'String',{'Parameter Sweep of GPe/STN Coupling and Resultant Local Beta Power'},...
%     'LineStyle','none',...
%     'FontSize',14,...
%     'FontWeight','bold',...
%     'FitBoxToText','off');




