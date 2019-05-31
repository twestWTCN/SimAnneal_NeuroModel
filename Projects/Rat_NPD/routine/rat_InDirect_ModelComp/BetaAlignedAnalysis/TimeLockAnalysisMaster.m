function TimeLockAnalysisMaster(R,BB,condsel,AS)
addpath(R.BBA_path)
TL.struccmap = BB.struccmap;
ip = 0;
figure
for cond = condsel
    ip = ip + 1;
    TL.periodT = [-800 800];
    TL = defineBurstTimeLockEpoch(BB,TL,cond);
    figure(3)
    subplot(1,3,ip)
    a = genBoxPlot(TL.onsetT{cond}',BB.struccmap,'horizontal');
    b= gca;
    xlim([-400 400])
    b.YTickLabel = R.chsim_name;
    ylabel('Region')
    xlabel('Burst Onset Time (ms)')
    grid on
    title(R.condname{cond})
    figure(4)
    subplot(2,3,ip)
    plotTLTimeEvolutions(TL,cond,'BP')
    title(R.condname{cond})
    xlabel('Burst Onset Time (ms)')
    xlim([TL.periodT])
    subplot(2,3,ip+3)
    plotTLTimeEvolutions(TL,cond,'amp',2.5)
    title(R.condname{cond})
    xlabel('Burst Onset Time (ms)')
    ylim([-10 0.5])
    xlim([TL.periodT])
%     figure(5)
%    plotTLPhaseSlipProb(TL,cond) 
%    hold on; 
end
set(gcf,'Position',[364         656        1445         322])
ip = 0;

% plot(TL.epochT,

for cond = condsel
    if isequal(R.condname{cond},'Empirical')
        T = BB.TEmp;
        TSw = BB.TSwEmp;
    else
        T = BB.T;
        TSw = BB.TSw;
    end
    ip = ip + 1;
    a(ip) = subplot(1,3,ip)
    for struc = 1:5
        if struc<5
            plot(T(1,1000:end-1000),1*zscore(BB.BPTime{cond}(struc,1000:end-1000)) - struc*10,'color',BB.struccmap(struc,:))
            %             plot(BB.T(1,1000:end-1000),0.1*(BB.BPTime{cond}(struc,1000:end-1000)) - struc*3,'color',cmap(struc,:))
        else
            plot(TSw(1,1000:end-1000),1*zscore(BB.PLV{cond}(1,1000:end-1000)) - struc*10,'color','b')
            %              plot(BB.TSw(1,1000:end-1000),0.1*(BB.PLV{cond}(1,1000:end-1000)) - struc*3,'color','b')
        end
        hold on
        xlim([120 122])
    end
end
linkaxes(a,'y')
set(F(5),'Position',[374.0000  88.0000   894.5000  224.0000])


F(6) = figure;
ip = 0;
for cond = condsel
    ip = ip + 1;
    subplot(1,3,ip)
    for struc = 1:4
        plot(R.frqz,AS{cond}(struc,:),'color',BB.struccmap(struc,:),'LineWidth',2); hold on
    end
    xlabel('Frequency (Hz)'); ylabel('Amplitude')
    title(R.condname{cond}); grid on
    ylim([0 35]);
end
set(F(6),'Position',[374.0000  88.0000   894.5000  224.0000])

close all
figure
a(1) = subplot(5,1,1)
X(1,:) = wrapToPi(BB.PhiTime{2}(1,:));
X(2,:) = wrapToPi(BB.PhiTime{2}(2,:));
X(3,:) = wrapToPi(BB.PhiTime{2}(4,:));
plot(T,X)

a(2) = subplot(5,1,2)
XR(1,:) = wrapToPi(BB.PhiTime{2}(1,:)-BB.PhiTime{1}(2,:));
XR(2,:) = wrapToPi(BB.PhiTime{2}(1,:)-BB.PhiTime{1}(4,:));
XR(3,:) = wrapToPi(BB.PhiTime{2}(2,:)-BB.PhiTime{1}(4,:));
plot(T,XR)
legend({'M1/STR','M1/STN','STR/STN'})

a(3) = subplot(5,1,3)
XR(1,:) = (BB.PhiTime{2}(1,:)-BB.PhiTime{1}(2,:));
XR(2,:) = (BB.PhiTime{2}(1,:)-BB.PhiTime{1}(4,:));
XR(3,:) = (BB.PhiTime{2}(2,:)-BB.PhiTime{1}(4,:));
XRD = diff(XR,[],2);
plot(T(2:end),XRD)

a(4) = subplot(5,1,4);
for struc = 1:5
    if struc<5
        plot(T(1,1000:end-1000),1*zscore(BB.BPTime{2}(struc,1000:end-1000)) - struc*10,'color',BB.struccmap(struc,:))
        %             plot(BB.T(1,1000:end-1000),0.1*(BB.BPTime{cond}(struc,1000:end-1000)) - struc*3,'color',cmap(struc,:))
    else
        plot(TSw(1,1000:end-1000),1*zscore(BB.PLV{2}(1,1000:end-1000)) - struc*10,'color','b')
        %              plot(BB.TSw(1,1000:end-1000),0.1*(BB.PLV{cond}(1,1000:end-1000)) - struc*3,'color','b')
    end
    hold on
    xlim([120 122])
end

a(5) = subplot(5,1,5)
plot(T, BB.data{2})

linkaxes(a,'x')