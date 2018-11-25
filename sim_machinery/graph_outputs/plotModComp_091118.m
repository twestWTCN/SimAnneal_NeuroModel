function plotModComp_091118(R,cmap)
% addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
% R.plot.confint = 'none';
if nargin<2
    cmap = brewermap(R.modcomp.modN,'Spectral');
    % cmap = linspecer(R.modcomp.modN);
end
%% First get probabilities so you can compute model space epsilon
for modID = 1:R.modcomp.modN
    dagname = sprintf([R.out.tag '_M%.0f'],modID);
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modeProbs_' R.out.tag '_NPD_' dagname '.mat'])
    A = varo; %i.e. permMod
    
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        r2bank{modID} = r2rep;
    end
end
r2bank = horzcat(r2bank{:});
R.modcomp.modEvi.epspop = prctile(r2bank,50); % threshold becomes median of model fits

p = 0;
mni = 0;
for modID = 1:R.modcomp.modN
    dagname = sprintf([R.out.tag '_M%.0f'],modID);
    load([R.rootn 'outputs\' R.out.tag '\NPD_' dagname '\modeProbs_' R.out.tag '_NPD_' dagname '.mat'])
    A = varo; %i.e. permMod
    
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2rep(isnan(r2rep) | isinf(r2rep)) = [];
        
        r2repSave{modID} = (r2rep);
        pd = fitdist(r2rep','Normal');
        x_values = -1:0.01:1;
        y = pdf(pd,x_values);
        figure(1)
        plot(x_values,y,'LineWidth',2)
        
        %         figure(1)
        %         histogram(r2rep,-1:0.1:1);
        hold on
        
        KL(modID) = sum(A.KL);
        DKL(modID) = sum(A.DKL);
        pmod(modID) = sum(r2rep>R.modcomp.modEvi.epspop)/ size(r2rep,2);
        
        h = figure(10);
        R.plot.cmap = cmap(modID,:);
        flag = 0;
        
        if ismember(modID,R.modcompplot.NPDsel)
            p = p +1;
            [hl(p), hp, dl, flag] = PlotFeatureConfInt_gen060818(R,A,h);
        end
        % hl(modID) = plot(1,1);
        if ~flag
            mni = mni +1;
            longlab{mni} = sprintf('Model %.f',modID);
        end
    else
        r2repSave{modID} = nan(size(r2rep));
    end
    
    shortlab{modID} = sprintf('M%.f',modID);
    %     pause(2)
    %     figure(10)
    %     h = findobj(gca,'Type','line');
    %     legend(h([size(h,1):-2:1 1]),longlab)
    
end
hl(end+1) = dl;
longlab{end+1} = 'Data';
figure(2)
subplot(3,1,1)
violin(r2repSave,'facecolor',cmap,'medc','k:','xlabel',shortlab)
hold on
plot([0 R.modcomp.modN+1],[R.modcomp.modEvi.epspop R.modcomp.modEvi.epspop],'k--') 
xlabel('Model'); ylabel('NMRSE'); grid on; ylim([-1.5 0.25])
a = gca; a.XTick = 1:R.modcomp.modN;
a.XTickLabel = shortlab;

figure(2)
subplot(3,1,1)
h = findobj(gca,'Type','line');
% legend(hl,{longlab{[R.modcompplot.NPDsel end]}})

subplot(3,1,2)
for i = 1:R.modcomp.modN
    b = bar(i,-log10(1-pmod(i))); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:R.modcomp.modN; grid on
a.XTickLabel = shortlab;
xlabel('Model'); ylabel('-log_{10} P(M|D)')

subplot(3,1,3)
for i = 1:R.modcomp.modN
    b = bar(i,KL(i)); hold on
    b.FaceColor = cmap(i,:);
end
a = gca; a.XTick = 1:R.modcomp.modN;
a.XTickLabel = shortlab;
grid on
xlabel('Model'); ylabel('KL Divergence')

% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')