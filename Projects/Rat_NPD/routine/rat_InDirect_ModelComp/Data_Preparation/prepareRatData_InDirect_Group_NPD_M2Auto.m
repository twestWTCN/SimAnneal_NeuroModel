function R = prepareRatData_InDirect_Group_NPD_M2Auto(R)
% prepareratdata_group(R.rootn,R.projectn);
load([R.filepathn '\NPD_paper_RatNPD_150618.mat']);
NPDmat = fyA;
load([R.filepathn '\nsPow_paper_RatNPD_150618.mat']);
nsPowmat = fyA;
load([R.filepathn '\frq_paper_RatNPD_150618.mat']);
F_data = fxA(:,1);
meannpd_data = [];
% condsel = [1 2];
R.condnames = {'OFF'};
R.Bcond = 1;
condsel = 2;
chsel = [1]'; % M1 STR GPe STN

for C =1:numel(R.condnames)
    %     X = squeeze(fyA(:,:,:,:,2,:)); {i,j,dirc,cond,sub}
    for i = 1:size(chsel,1)
        for j = 1:size(chsel,1)
            if i==j
                Ftmp = F_data;
                Pxy = abs(log10(mean(vertcat(nsPowmat{chsel(i),condsel(C),:}),1)));
                %                 Pxy = Pxy(Ftmp>=R.frqz(1) & Ftmp<=R.frqz(end));
                %                 Ftmp = Ftmp(Ftmp>=R.frqz(1) & Ftmp<=R.frqz(end));
                Pxy(Ftmp>48 & Ftmp<52) = [];
                Ftmp(Ftmp>48 & Ftmp<52) = [];
                [xCalc yCalc b Rsq] = linregress(log10(Ftmp),Pxy');
                Pxy = Pxy-yCalc';
                Pxy = interp1(Ftmp,Pxy,R.frqz);
%                 plot(R.frqz,Pxy); hold on
                f = fit(R.frqz',Pxy','gauss3');
                Pxy = f(R.frqz)';
%                 plot(R.frqz,Pxy); hold on
                Pxy = 10.^Pxy;
                Pxy = (Pxy-mean(Pxy))./std(Pxy);
                Pxy = Pxy - min(Pxy);
                %                 Pxy = Pxy.*tukeywin(length(Pxy),0.6)';
                
                %                 plot(R.frqz,10.^(Pxy));
                meannpd_data(C,i,j,1,:) = Pxy;
%                 close all
            else
                
                for k = 1:size(NPDmat,3)
                    Ftmp = F_data;
                    Cxy = mean(horzcat(NPDmat{chsel(i),chsel(j),k,condsel(C),:})',1);
                    Cxy = Cxy*1.5;
                    Cxy(Ftmp>48 & Ftmp<52) = [];
                    Ftmp(Ftmp>48 & Ftmp<52) = [];
                    Cxy = interp1(Ftmp,Cxy,R.frqz);
                    %                     Cxy = Cxy.*tukeywin(length(Cxy),0.6)'; %%NPD_sim_n(i,j,1,:)
%                     plot(R.frqz,Cxy); hold on
                    f = fit(R.frqz',Cxy','gauss3');
                    Cxy = f(R.frqz)';
%                     plot(R.frqz,Cxy)
                    meannpd_data(C,i,j,k,:) = Cxy;
%                     close all
                end
                
            end
        end
    end
end

% Set data as working version
R.data.feat_emp = meannpd_data;
% squeeze(meannpd_data(1,1,1,1,:))
R.data.feat_xscale = R.frqz;

% Plot CSD
if strcmp('CSD',R.data.datatype)
    csdplotter_220517({meannpd_data},[],F_data,R)
elseif strcmp('NPD',R.data.datatype)
    npdplotter_110717({meannpd_data},[],R.frqz,R,[],[])
end