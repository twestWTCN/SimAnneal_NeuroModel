function [] = optProgPlot(tm,tbr2,p,R)
    subplot(2,1,1)
    plot(-tm,tbr2)
    ylim([-1 1])
    subplot(2,1,2)
    a = fieldnames(p);
    for i = 1:numel(a)
%         if numel(eval(['p.' a{i}]))>100 &&
        if strmatch(a{i}(end),'_s')
            p = rmfield(p,a{i});
        end
    end
    par = full(spm_vec(p));
    par(par<-20) = 0;
    bar(par)
%     ylim([1.2*min(par) 1.2*max(par)]);
    xlim([0 length(par)])
    ylim([R.SimAn.pOptBound])
    xlabel('parameter')
    ylabel('Posterior')