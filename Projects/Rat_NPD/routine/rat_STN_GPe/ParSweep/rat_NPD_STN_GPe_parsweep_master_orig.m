clear ; close all

%  simAnnealAddPaths()
rng(23123)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
modID = 1;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
R.analysis.modEvi.eps = prctile(parBank(end,:),75); % Get models above 75th percentile
R.analysis.modEvi.N = 1000;
optP = getOptParMean(m,p,R,parBank);

%% Parameter Sweep across STN/GPe Subcircuit Connections
N = 8;
parsweep.N = N;
parsweep.R = '.A{1}(1,2)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Q = '.A{2}(2,1)';
parsweep.Qname = 'GPe-|STN Connections Strength';
parsweep.Qlist = linspace(-2,3,N);

parsweep.InvertXY = [optP.A{1}(1,2) optP.A{2}(2,1)];

% Conduct sweep
R.obs.gainmeth = {'unitvar'};
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R = simannealsetup_NPD_STN_GPe;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
pathstr = [R.rootn 'outputs\' R.out.tag '\' R.out.dag];
mkdir(pathstr); 
save([pathstr '\'  R.out.dag '_A_parsweep'],'parsweep');
load([pathstr '\'  R.out.dag '_A_parsweep'],'parsweep');

% Plot heatmap
figure(4)
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
    elseif i == 2
        parsweep.plot.feat = log10(squeeze(parsweep.frqPowBank(2,:,:)));
    end
    subplot(2,2,i)
    parSweepPlot(R,parsweep)
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    if i==1
        caxis([5 21])
    elseif i==2
        caxis([-1 0])
    end
end

%% Parameter Sweep Time CoOnstants
% Reload Data
R = simannealsetup_NPD_STN_GPe;
modID = 1;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
R.analysis.modEvi.eps = -0.5;
R.analysis.modEvi.N = 1000;
optP = getOptParMean(m,p,R,parBank);

N = 128;
parsweep.N = N;
parsweep.R = '.int{1}.T';
parsweep.Rname = 'GPe Time Constant';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Q = '.int{2}.T';
parsweep.Qname = 'STN Time Constant';
parsweep.Qlist = linspace(-2,3,N);
parsweep.InvertXY = [optP.int{1}.T optP.int{2}.T];

% Conduct sweep
R.obs.gainmeth = {'unitvar'};
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
pathstr = [R.rootn 'outputs\' R.out.tag '\' R.out.dag];
mkdir(pathstr); 
save([pathstr '\'  R.out.dag '_STN_GPe_parsweep'],'parsweep');
load([pathstr '\'  R.out.dag '_STN_GPe_parsweep'],'parsweep');

% Plot heatmap
figure(4)
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
    elseif i == 2
        parsweep.plot.feat = log10(squeeze(parsweep.frqPowBank(2,:,:)));
    end
    subplot(2,2,i+2)
    parSweepPlot(R,parsweep)
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    if i==1
        caxis([5 21])
    elseif i==2
        caxis([-1 0.1])
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




