function [R] = getFlaviesData(R)
load('GA_data_model.mat')
% Convert to rads!
Data.Phase_Amp(:,1) = deg2rad(Data.Phase_Amp(:,1));
Data.Phase_AmpCV(:,1) = deg2rad(Data.Phase_AmpCV(:,1));
Data.Phase_Lat(:,1) = deg2rad(Data.Phase_Lat(:,1));
Data.Phase_LatCV(:,1) = deg2rad(Data.Phase_LatCV(:,1));

% conver %CV to CV
Data.Phase_AmpCV(:,2) = Data.Phase_AmpCV(:,2)./100;
Data.Phase_LatCV(:,2) = Data.Phase_LatCV(:,2)./100;

feat_emp{1} = Data.Phase_Amp';
feat_emp{2} = Data.Phase_AmpCV';
feat_emp{3} = Data.Phase_Lat';
feat_emp{4} = Data.Phase_LatCV';

R.trans.amprange = -3:0.25:3.5;
x = Data.AmpvsLat(:,1);
x = (x-mean(x))./std(x);

[binStats_ampLag] = binDatabyRange(x,R.trans.amprange,Data.AmpvsLat(:,2).*1000);
binStats_ampLag(end,:) = [];
y = binStats_ampLag(:,3)'; % Binned Lags by Amp
y(isnan(y)) = 0; % Fill nans with mean
feat_emp{5} = [binEdge2Mid(R.trans.amprange); y];

R.trans.pdfSupport = linspace(-4,4,32);
y =Data.AmpvsLat(:,1);
y = (y-mean(y))./std(y);
x =Data.AmpvsLat(:,2).*1000;
x1 = x(y>1.2);
[a] = histcounts(x1,R.trans.pdfSupport,'Normalization','probability');
feat_emp{6} = [R.trans.pdfSupport(1:end-1);a];
plot(feat_emp{6}(1,:),feat_emp{6}(2,:),'LineWidth',1,'color',[0 0 1],'LineStyle','--')

x2 = x(y<=1.2);
[a] = histcounts(x2,R.trans.pdfSupport,'Normalization','probability');
feat_emp{7} = [R.trans.pdfSupport(1:end-1);a];
plot(feat_emp{7}(1,:),feat_emp{7}(2,:),'LineWidth',1,'color',[0 0 1],'LineStyle','--')


% feat_emp{6} = [1000.*Data.Lat_Distrib_NoOpt_Phase(:,1)';Data.Lat_Distrib_NoOpt_Phase(:,2)'];
% feat_emp{7} = [TMS_phase; MEP_amp; MEP_onset];
R.data.feat_emp = feat_emp;
% R.trans.pdfSupport = Data.Lat_Distrib_NoOpt_Phase(:,1).*1000';