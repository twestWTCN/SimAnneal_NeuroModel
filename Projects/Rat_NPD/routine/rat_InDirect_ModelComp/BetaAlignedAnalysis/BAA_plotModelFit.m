function [R] = BAA_plotModelFit(Rorg,modID,simtime)
Rorg.obs.gainmeth = {'unitvar','boring'};
 Rorg.obs.trans.norm = 1; Rorg.obs.trans.gauss = 1;
 Rorg.obs.obsstates = [1:6];
 Rorg.chloc_name = Rorg.chsim_name;
[R,m,permMod,xsimMod] = getSimModelData_v2(Rorg,modID,simtime,1);


close all; 
X = xsimMod{1}{1}{1};
for i = 1:6
    figure(1)
    plot(R.IntP.tvec_obs,X(i,:))
    hold on
%     figure(2+i)
%  pspectrum(X(i,:),1/R.IntP.dt,'spectrogram','TimeResolution',0.5,...
% 'OverlapPercent',99,'MinTHreshold',-20);
% ylim([0 0.098])
    
end
xlim([30 40])
set(gcf,'Position',[50         550        1686         436])
R = prepareRatData_NoGauss_Group_NPD(R);
figure
plotABCSpectraOnly(R.data.feat_xscale,R.data.feat_emp,permMod{1}.feat_rep{1})
figure
ab = npdplotter_110717({R.data.feat_emp},{permMod{1}.feat_rep{1}},R.data.feat_xscale,R,[],[])

for i = 1:36
    if ~any(intersect([1 8 15 22 29 36],i))
        subplot(6,6,i)
        ylim([0 1])
    end
end


% A = X(6,:);
% 
% STN_slow = ft_preproc_bandpassfilter(X(4,:),2000,[0.5 4],[],'fir');
% 
% slow_list = linspace(2.5,6,32);
% fast_list = linspace(8,34,32);
% Rz = [];
% for i = 1:size(slow_list,2)
%     for j = 1:size(fast_list,2)
%         A = ft_preproc_bandpassfilter(X(4,:),2000,[slow_list(i)-2 slow_list(i)+2],[],'fir');
%         B = abs(hilbert(ft_preproc_bandpassfilter(X(4,:),2000,[fast_list(j)-2 fast_list(j)+2],[],'fir')));
%         
%        LR = corrcoef(A,B);
%        Rz(i,j) = LR(2);
%        if i==4
%            [r,lags] = xcorr(A,B,5000);
%            xcorstore(:,j) = r;
%        end
%         disp([i j])
%     end
% end
% subplot(1,3,1)
% imagesc(slow_list,fast_list,Rz')
% set(gca,'YDir','normal');
% subplot(1,2,2)
% imagesc(lags,fast_list,xcorstore')
% set(gca,'YDir','normal'); xlim([-1000 1000])
% % subplot(1,3,3)
% % imagesc(lags,fast_list,lagstore_j')
% % set(gca,'YDir','normal'); xlim([-1000 1000])
% 
% 
% 
% 
