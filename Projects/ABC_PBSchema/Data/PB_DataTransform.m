function [dum,feat,wflag] = PB_DataTransform(R,xsims)
wflag = 0; dum = [];
% R.trans.pirange = linspace(-pi,pi,8);
% feat_emp{1} = [binEdge2Mid(xrange);binStats_phiAmp(:,3)'];
% feat_emp{2} = [binEdge2Mid(xrange);binStats_phiLag(:,3)'];
% R.trans.lagrange = 0.5:0.75:35;
% feat_emp{3} = [binEdge2Mid(xrange);binStats_ampLag(:,3)'];
% feat_emp{4} = [lagX;lagPDF];

[binStats_phizAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,xsims.zMEP_amp);
binStats_phizAmp(end,:) = [];
[y,K] = recentreDist(binStats_phizAmp(:,3));
binStatsRS_phizAmp = y;

[binStats_phizLag] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,xsims.zMEP_onset);
binStats_phizLag(end,:) = [];
binStatsRS_phizLag = circshift(binStats_phizLag(:,3),K,1);

% Take non-Zscored for CoV (no negative values!)

[binStats_phiAmp] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,xsims.MEP_amp);
binStats_phiAmp(end,:) = [];
binStatsRS_phiAmpCoV = circshift(binStats_phiAmp(:,4),K,1).*10;

[binStats_phiLag] = binDatabyRange(xsims.TMS_phase,R.trans.pirange,xsims.MEP_onset);
binStats_phiLag(end,:) = [];
binStatsRS_phiLagCoV = circshift(binStats_phiLag(:,4),K,1).*10;

% Amp by mean corrected lag
[binStats_ampLag] = binDatabyRange(xsims.zMEP_amp,R.trans.amprange,xsims.MEP_mconset);
binStats_ampLag(end,:) = [];

% [lagX,lagPDF] = fitBimodal(xsims.MEP_onset,R.trans.pdfSupport,[.5 -1 -1 1 1],[0 -5 -5 0.001 0.001],[1 200 200 25 25]);
y =xsims.zMEP_amp(1,:);
x =xsims.MEP_mconset(1,:);
x1 = x(y>1.2);
lagPDF = histcounts(x1,R.trans.pdfSupport,'Normalization','probability');

xrange = R.trans.pirange;
feat{1} = [binEdge2Mid(xrange);binStatsRS_phizAmp']; % Z Amplitude to phase
feat{2} = [binEdge2Mid(xrange);binStatsRS_phiAmpCoV']; % CoV Amplitude to phase
feat{3} = [binEdge2Mid(xrange);binStatsRS_phizLag'];    % Z Lag to phase
feat{4} = [binEdge2Mid(xrange);binStatsRS_phiLagCoV']; % CoV Lag to phase
xrange = R.trans.amprange;

y = binStats_ampLag(:,3)';
y(isnan(y)) = 0; % Fill nans with mean
feat{5} = [binEdge2Mid(xrange);y]; % Lag to Amp
feat{6} = [R.trans.pdfSupport(1:end-1); lagPDF];            % Lag Bimodal PDF

y =xsims.zMEP_amp(1,:);
x =xsims.MEP_mconset(1,:);
x2 = x(y<=1.2);
lagPDF = histcounts(x2,R.trans.pdfSupport,'Normalization','probability');
feat{7} = [R.trans.pdfSupport(1:end-1); lagPDF];            % Lag Bimodal PDF

% feat{7} = [xsims.TMS_phase; xsims.MEP_amp; xsims.MEP_onset];
