clear ; close all

%  simAnnealAddPaths()
rng(23123)

%% Set Routine Pars
R = simannealsetup_NPD_STN_GPe;
modID = 1;
R.out.dag = sprintf(['NPD_' R.out.tag '_M%.0f'],modID);
[R,m,p,parBank] = loadABCData(R);
R.analysis.modEvi.eps = -0.5;
R.analysis.modEvi.N = 1000;
optP = getOptParMean(m,p,R,parBank);

%% Parameter Sweep across STN/GPe Subcircuit Connections
N = 48;
parsweep.N = N;
parsweep.R = '.A{1}(3,4)';
parsweep.Rname = 'STN-> GPe Connections Strength';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Q = '.A{2}(4,3)';
parsweep.Qname = 'GPe-|STN Connections Strength';
parsweep.Qlist = linspace(-2,3,N);

% Conduct sweep
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R = simannealsetup_NPD_STN_GPe;
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STN_GPe_TC_parsweep'],'parsweep');

% Plot heatmap
figure(4)
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
    elseif i == 2
        parsweep.plot.feat = squeeze(parsweep.betaPowBank(2,:,:));
    end
    subplot(2,2,i+1)
    parSweepPlot(R,parsweep)
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    for i=~1
        % caxis([0 0.1])
    end
end

%% Parameter Sweep Time CoOnstants
N = 32;
parsweep.N = N;
parsweep.R = '.int{1}.T';
parsweep.Rname = 'GPe Time Constant';
parsweep.Rlist = linspace(-2,3,N);
parsweep.Q = '.int{2}.T';
parsweep.Qname = 'STN Time Constant';
parsweep.Qlist = linspace(-2,3,N);

% Conduct sweep
parsweep = modelBetaParSweep(m,optP,parsweep,R);
R = simannealsetup_NPD_STN_GPe;
pathstr = [R.rootn 'analysis\parsweeps\'];
mkdir(pathstr);
save([pathstr '\A_STN_GPe_TC_parsweep'],'parsweep');

% Plot heatmap
figure(4)
for i = 1:numel(R.chsim_name)
    if i == 1
        parsweep.plot.feat = squeeze(parsweep.maxfrqBank(2,:,:));
    elseif i == 2
        parsweep.plot.feat = squeeze(parsweep.betaPowBank(2,:,:));
    end
    subplot(2,2,i+2)
    parSweepPlot(R,parsweep)
    xlabel(parsweep.Rname,'FontSize',10)
    ylabel(parsweep.Qname,'FontSize',10)
    axis square
    title([R.chsim_name{2}],'FontSize',10,'fontweight','bold')
    c = colorbar;
    for i=~1
        % caxis([0 0.1])
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




