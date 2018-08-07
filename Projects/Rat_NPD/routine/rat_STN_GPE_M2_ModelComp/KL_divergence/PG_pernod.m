function plotModComp(R,permMod)
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\violin')
R.plot.confint = 'none';
for mnum = 1:7
    A = permMod{mnum};
    
    if ~isempty(A)
        r2rep = [A.r2rep{:}];
        r2repSave{mnum} = r2rep;
        pd = fitdist(r2rep','Normal');
        x_values = -1:0.01:1;
        y = pdf(pd,x_values);
        figure(1)
        plot(x_values,y,'LineWidth',2)
        
        %         figure(1)
        %         histogram(r2rep,-1:0.1:1);
        hold on
        if mnum ~= 7
            KL(mnum) = sum(A.KLrep{1});
            DKL(mnum) = sum(A.DKLrep{1});
        else
            KL(mnum) = sum(A.KL);
            DKL(mnum) = sum(A.DKL);
        end
        pmod(mnum) = sum(r2rep>-0.25)/1000;
        h = figure(10);
        R.plot.cmap = cmap(mnum,:);
        [hl, hp] = PlotFeatureConfInt_gen060818(R,permMod{mnum},h)
    else
        r2repSave{mnum} = nan(size(r2rep));
    end
end
figure(30)
violin(r2repSave,'facecolor',cmap,'medc','k:')
a = gca; a.XTickLabel = {'M1','M2','M3','M4','M5','M6','M7'};

figure(10)
h = findobj(gca,'Type','line');
legend(h(12:-2:1),'Model 1','Model 2','Model 4','Model 5','Model 6','Model 7')
figure(2)
subplot(2,1,1)
bar(pmod)
xlabel('Model'); ylabel('P(M|D)')
subplot(2,1,2)
bar(KL)
xlabel('Model'); ylabel('KL Divergence')
% subplot(3,1,3)
% bar(DKL)
% xlabel('Model'); ylabel('Joint KL Divergence')