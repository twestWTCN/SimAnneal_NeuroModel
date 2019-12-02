function [R] = BAA_sim_lesionExp_cortical_rhythms(R,simtime,fresh)
close all
% Comopute simulations by sweeping across data
[R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,10,simtime);
p = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);
R.Bcond = -1;
cmap = brewermap(128,'RdBu');
ckeypow = linspace(-100,75,128);
ckeyfrq = linspace(-10,10,128);

% PRC.condcmap = cmap([1 4 8 16 4 18],:);

% G(:,1)  mp -> mp (-ve self)  4
% G(:,2)  mp -> sp (+ve rec )  4
% G(:,3)  ii -> mp (-ve rec )  4
% G(:,4)  ii -> ii (-ve self)  4
% G(:,5)  mp -> ii (+ve rec )  4
% G(:,6)  dp -> ii (+ve rec )  2
% G(:,7)  sp -> sp (-ve self)  4
% G(:,8)  sp -> mp (+ve rec )  4
% G(:,9)  ii -> dp (-ve rec )  2
% G(:,10) dp -> dp (-ve self)  1
% G(:,11) sp -> dp (+ve rec)  2
% G(:,12) ii -> sp (-ve rec)  4
% G(:,13) sp -> ii (+ve rec)  4
% G(:,14) dp -> sp (+ve rec)  2
%---------------------------------

condname = {'base','No CTX -> Thal','No Thal -> CTX','No GPi->Thal',...
    'mp -> mp (-ve self)','mp -> sp (+ve rec )','ii -> mp (-ve rec )',...
    'ii -> ii (-ve self)','mp -> ii (+ve rec )','dp -> ii (+ve rec )',...
    'sp -> sp (-ve self)','sp -> mp (+ve rec )','ii -> dp (-ve rec )',...
    'dp -> dp (-ve self)','sp -> dp (+ve rec)','ii -> sp (-ve rec)',...
    'sp -> ii (+ve rec)','dp -> sp (+ve rec)'};

MP.m = m;
MP.p = p;
% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
R.obs.trans.norm = 0;
% R.obs.outstates = find([MP.m.outstates{:}]);
% R.obs.obsstates = 1:4;
Pbase = p;
% powIJ_B = {}; peakIJ_B = {}; freqIJ_B{};
if fresh == 1
    for i = 5:size(condname,2)
        Pbase_i = Pbase;
        %
        if i==2
            Pbase_i.A{1}(6,1) = -128; % Sever the thalamic input
        elseif i==3
            Pbase_i.A{1}(1,6) = -128; % Sever the thalamic input
        elseif i==4
            Pbase_i.A{2}(6,5) = -128; % Sever the thalamic input       
        elseif i>4
            Pbase_i.A{1}(1,6) = -128; % Sever the thalamic input
            Pbase_i.int{1}.G(i-4) = -128;
        end
        [r2,~,feat_sim] = computeSimData(R,m,uc,Pbase_i,0);
        [powIJ_B(i),peakIJ_B(i),freqIJ_B(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,3,:)),[14 30]);
        [powIJ_B1(i),peakIJ_B1(i),freqIJ_B1(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,3,:)),[14 21]);
        [powIJ_B2(i),peakIJ_B2(i),freqIJ_B2(i)] = findSpectralStats(R.frqz,squeeze(feat_sim(1,1,1,3,:)),[21 30]);
        fitIJ(i) = r2;
        disp([i])
        figure(1)
        plot(R.frqz,squeeze(feat_sim(1,1,1,3,:)));
        hold on
        
    end
    save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion_CTX'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
else
    load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion_CTX'],'powIJ_B','peakIJ_B','freqIJ_B',...
        'powIJ_B1','peakIJ_B1','freqIJ_B1',...
        'powIJ_B2','peakIJ_B2','freqIJ_B2','condname')
end
% Plot heatmaps
for band = 1:3
    figure
    titlist = {'14-30 Hz','14-21 Hz','21-30 Hz'}
    if band == 1
        X = (powIJ_B-powIJ_B(1))./(powIJ_B(1)).*100;
        Y = (freqIJ_B-freqIJ_B(1));
    elseif band == 2
        X = (powIJ_B1-powIJ_B1(1))./(powIJ_B1(1)).*100;
        Y = (freqIJ_B1-freqIJ_B1(1));
    elseif band == 3
        X = (powIJ_B2-powIJ_B2(1))./(powIJ_B2(1)).*100;
        Y = (freqIJ_B2-freqIJ_B2(1));
    end
    bcmap = brewermap(128,'RdGy');
    subplot(1,2,1)
    colormap(bcmap)
    igsc = imagesc(X)
    set(gca,'YDir','normal');
    cb = colorbar;
    caxis([-100 100])
    %         title(titn)
    
    % Plot as Barplots
    subplot(2,1,1)
    b = bar(X);
    [dum keyind] = min(abs(X'-ckeypow),[],2);
    %     b.CData = cmap(keyind,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('% Power Change')
    %     xlim([1.5 11.5])
    subplot(2,1,2)
    b = bar(Y);
    [dum keyind] = min(abs(Y'-ckeyfrq),[],2);
    %     b.CData = cmap(keyinBd,:);
    b.EdgeAlpha = 0; b.FaceColor = 'flat';
    a = gca; grid on; box off
    a = gca;
    a.XTickLabel = condname; a.XTickLabelRotation = 45; ylabel('Change in Peak Frequency')
    %     xlim([1.5 11.5])
    % Make Tables
    tabdata(:,band) = X;
    
end
set(gcf,'Position',[448   121   559   592])
%
% Z = (powIJ_B-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the real (possibly nonlinear effect)
% Z(getDiagInd(Z)) = zeros(1,size(Z,1));
% X = powIJ_B;
% Y = X;
% for i = 1:size(X,1)
%     for j= 1:size(X,2)
%         if i~=j
%             Y(i,j) = ((X(i,i)+X(j,j))-powIJ_B(1,1))./(powIJ_B(1,1)).*100; % This is the linear effect
%         end
%     end
% end
%
% figure
% colormap(flipud(cmap))
% Q = abs(Y)-abs(Z)
% Z(isnan(Z)) = 0;
% igsc = imagesc(Q)
% set(gca,'YDir','normal');
% caxis([-125 125]);
% a = gca; a.XTick = 1:size(X,1); a.XTickLabel = condname; a.XTickLabelRotation = 45;
% a.YTickLabel = condname; a.YTickLabelRotation = 0;
% cb = colorbar; axis equal; xlim([1.5 11.5]);ylim([0.5 11.5])
%
% set(gcf,'Position',[1065          84         692         592])
%
%
% % subplot(1,2,2)
% % colormap(cmap)
% % igsc = imagesc(fitIJ)
% % set(gca,'YDir','normal');
% % cb = colorbar;
% %
% save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\BAA_lesion'],'powIJ_B','powIJ_B1','powIJ_B2','condname')

% plotSimulatedPRC(PRC,[2 1 3])
% set(gcf,'Position',[711         604        1081         374])