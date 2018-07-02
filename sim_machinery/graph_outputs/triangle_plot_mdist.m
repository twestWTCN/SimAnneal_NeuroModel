function triangle_plot_mdist(R,p,m,xf)
% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m,1); % in structure form
pIndMap = spm_vec(pInd); % in flat form

r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,5000);

labellist = {'M1 -> STR','M1 -> STN','STN -> GPe','STN -> GPi','Thal -> M1'};
set(gcf,'Position',[667         67        1024         1024])
%% INDEXING ISNT RIGHT FOR r THIS IS THE OPTIMIZED
A1 =  pInd.A{1}; %reshape(pInd.A{1},6,6);
plist(1) = A1(1); % M1 to STR
plist(2) = A1(2); % M1 to STN
plist(3) = A1(3); % STN to GPe
plist(4) = A1(4); % STN to GPi
plist(5) = A1(5); % THAL to M1
plist(6) = A1(6);
plist(7) = A1(7);
plist(8) = A1(8);

% A2 =  reshape(pInd.A{2},6,6);
% plist(6) = A2(3,2);
% plist(7) = A2(5,2);
% plist(8) = A2(4,3);
cnt = 0;

for j = 1:5
    for i = 1:5
        cnt = cnt+1;
        if i<=j
            I = find(pIndMap==plist(i));
            J = find(pIndMap==plist(j));
            
            u1 = r(:,I);
            v1 = r(:,J);
            x1 = ksdensity(xf(I,:),u1,'function','icdf');
            y1 = ksdensity(xf(J,:),v1,'function','icdf');
            
            [n,c]  = hist3([x1, y1],[15 15]);
            % scatter(x1,y1);
            subplot(5,5,cnt)

            [dum a] = contour(c{1},c{2},n,10);
            xlim(R.SimAn.pOptBound.*0.075)
            ylim(R.SimAn.pOptBound.*0.075)
            a = gca;
            a.FontSize = 14;

            
            if i==1
                ylabel(labellist{j});
            else
                a.YTickLabel = [];
            end
            if j==5
                xlabel(labellist{i});
            else
                a.XTickLabel = [];
            end
              grid on          
        end
    end
end

annotation(gcf,'textbox',...
    [0.27834375 0.927734375 0.4921640625 0.0576171875],...
    'String',{'Sample of Multivariate Posterior Distribution of Connectivity Parameters'},...
    'LineStyle','none',...
    'HorizontalAlignment','center',...
    'FontSize',18,...
    'FitBoxToText','off');