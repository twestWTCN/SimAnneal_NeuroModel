function BA_analysis(R,data,permMod)
ftdata = [];
ftdata.label = R.chsim_name;
ftdata.trial{1} = data{1};
ftdata.time{1} = R.IntP.tvec_obs;
ftdata.fsample = 1/R.IntP.dt;

cfg = [];
cfg.length = 1;
ftdata_tr = ft_redefinetrial(cfg,ftdata);

AS = squeeze(permMod.feat_rep{1}(1,:,1,1,:)); % get autospectra
banddef = [18 21]; % define band
frq_peak = getMaxBandFrq(AS,banddef,R.frqz);

% Compute Wavelet Amplitude
% cfg = [];
% cfg.method     = 'wavelet';
% cfg.width      = 10; % spectral bandwidth at frequency F: (F/width)*2
% cfg.gwidth     = 2;
% cfg.output     = 'pow';
% cfg.foi        = 14:0.5:35;
% cfg.toi        = ftdata.time{1};
% cfg.pad = 100;
% TFRwave = ft_freqanalysis(cfg, ftdata);
% fbwid = (median(banddef(1,:))/cfg.width)*2;
% disp(['Beta amp estimate at mid freq ' num2str(TFRwave.freq(fInd)) ' Hz +/- ' num2str(fbwid)])
%

% Optional Plots
% cfg = []
% cfg.channel = 'MMC';
% % cfg.zlim = [0 60]
% ft_singleplotTFR(cfg, TFRwave)


% Define TimeSeries
% % BB.Time = TFRwave.time;
BB.Time = ftdata.time{1};
BB.fsamp = 1/diff(BB.Time(end-1:end));
% [dum fInd] = min(abs(TFRwave.freq-frq_peak(1,4))); % Take middle!
% fbwid = (median(banddef(1,:))/cfg.width)*2;
% disp(['Beta amp estimate at mid freq ' num2str(TFRwave.freq(fInd)) ' Hz +/- ' num2str(fbwid)])


% % BB.Amp = squeeze(TFRwave.powspctrm(:,fInd,:));
% % BB.epsAmp = prctile(BB.Amp(4,:),75,2);

% interpolate original series to new sample
BB.RawInt = [];  BB.Amp = [];
for i = 1:size(ftdata.trial{1},1)
    bdata = ft_preproc_bandpassfilter(ftdata.trial{1}(i,:),ftdata.fsample,[frq_peak(1,4)-1 frq_peak(1,4)+1],[],'fir');
    BB.RawInt(i,:) = bdata;
    % %     BB.RawInt(i,:) = interp1(ftdata.time{1},bdata,BB.Time);
    BB.Amp(i,:) = abs(hilbert(BB.RawInt(i,:)));
end

BB.epsAmp = prctile(BB.Amp(4,:),75,2);

X = BB.Amp(4,:);
Xcd = X>BB.epsAmp;
Xcd = double(Xcd);

period = (3/banddef(1,1))*BB.fsamp;
consecSegs = SplitVec(find(Xcd(1,:)),'consecutive');

% Work first with lengths
segL = cellfun('length',consecSegs);
segInds = find(segL>(period)); % segs exceeding min length
segL_t = (segL/BB.fsamp)*1000; % Segment lengths in ms
segL_t(setdiff(1:length(segL),segInds)) = [];

% Segment Time indexes
segT_ind = [];
for ci = 1:numel(segInds); segT_ind{ci} = (consecSegs{segInds(ci)}([1 end])/BB.fsamp); end

% Now do Amplitudes
segA = [];
for ci = 1:numel(segInds); segA(ci) = nanmean(X(consecSegs{segInds(ci)})); end

cmap = linspecer(4); %brewermap(4,'Set1');
cmap = cmap(end:-1:1,:);
onsetT = [-500 500];
BEpoch = []; REpoch = [];
for i = 1:numel(segInds)
    epochdef = [consecSegs{segInds(i)}(1)+ floor((onsetT(1)/1e3)*BB.fsamp):consecSegs{segInds(i)}(1) + floor((onsetT(2)/1e3)*BB.fsamp)];
    if epochdef(end)<size(BB.Amp,2) && epochdef(1) > 0
        BEpoch(:,:,i) = BB.Amp(:,epochdef);
        REpoch(:,:,i) = BB.RawInt(:,epochdef);
    end
end
TEpoch = linspace(onsetT(1),onsetT(2),size(epochdef,2));
meanEnv = squeeze(mean(BEpoch,3)); stdEnv = squeeze(std(BEpoch,[],3)); %./sqrt(size(BEpoch,3)));
meanRaw = squeeze(mean(REpoch,3));
[x ind] = max(meanEnv,[],2)
for i = 1:size(BEpoch,1)
    XY = 3.5*meanEnv(i,:);
    XYp = XY-mean(XY);
    XY = 10 + ((XYp)-(i*2));
    lhan((2*i)-1) = plot(TEpoch,XY','color',cmap(i,:),'LineWidth',2); hold on
    [lp hp] = boundedline(TEpoch,XY',stdEnv(i,:)')
    lp.Color = cmap(i,:);
    hp.FaceColor = cmap(i,:); hp.FaceAlpha = 0.5;
    XY = 3*(10.^XYp).*meanRaw(i,:);
    XY = XY - mean(XY);
    XY = 10 + ((XY)-(i*2));
    plot(TEpoch,XY','color',cmap(i,:))
    lhan(2*i) = plot(repmat(TEpoch(ind(i)),2,1),[1.5 500],'LineStyle','--','color',cmap(i,:),'LineWidth',1);
    splay = TEpoch(ind(i))+(10*(TEpoch(ind(i))-median(TEpoch(ind))));
    plot([splay TEpoch(ind(i))],[0.5 1.5],'LineStyle','--','color',cmap(i,:),'LineWidth',1)
    text(splay,0.4,sprintf('%0.1f ms',TEpoch(ind(i))),'color',cmap(i,:))
    
    legn{(2*i)-1} = R.chsim_name{i};
    legn{2*i} = [R.chsim_name{i} ' max'];
end
% legend(lhan(1:2:8),legn(1:2:8))

xlabel('STN Beta Onset Time (ms)'); ylabel('Average Amplitude (a.u.)'); ylim([0 10]); xlim(onsetT)
set(gca,'YTickLabel',[]);
title('Sequential Beta Onset in Model 5')
set(gcf,'Position',[680          85         419        1013])