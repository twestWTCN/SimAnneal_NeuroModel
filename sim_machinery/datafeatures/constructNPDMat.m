function [F meannpd] = constructNPDMat(data,chloc_name,chlist,fsamp,N,R)
if isempty(N)
    N = floor(fsamp/R.obs.csd.df);
end
% N = R.obs.csd.pow2;
F_scale = R.frqz;
% Construct NPD matrix from data - take mean across channel replicates
for chloc = 1:size(chloc_name,2)
    chinds{chloc} = strmatch(chloc_name{chloc},chlist);
end
for chI = 1:size(chloc_name,2)
    for chJ = 1:size(chloc_name,2)
        for p = 1:size(chinds{chI},1)
            chindsP = chinds{chI};
            for r = 1:size(chinds{chJ},1)
                chindsR = chinds{chJ};
                if chI == chJ
                    [Pxy,F] = pwelch(data(chindsP(p),:),hanning(2^N),[],2^(N),fsamp);
                    Pxy = -abs(log10(Pxy));
                    Pxy(F>48 & F<52) = [];
                    F(F>48 & F<52) = [];
                    Pxy(F==0) = [];
                    F(F==0) = [];
                    [xCalc yCalc b Rsq] = linregress(log10(F),Pxy);
                    Pxy = Pxy-yCalc;
                    if nargin>5
                        Pxy = interp1(F,Pxy,F_scale);
                    else
                        Pxy =  Pxy(F>4);
                    end
                    Pxy = (Pxy-mean(Pxy))./std(Pxy);
                    Pxy = Pxy - min(Pxy);
                    Pxy = Pxy; %.*tukeywin(length(Pxy),0.25)';
                    xcsd(p,r,1:3,:) = repmat(Pxy,3,1);
                else
                    [f13,~,~]=sp2a2_R2(data(chindsP(p),:)',data(chindsR(r),:)',fsamp,N-1);
                    F = f13(:,1);
                    zl = [10 11 12];
                    for z = 1:3
                        %                     [Pxy,F] = cpsd(data(chindsP(p),:),data(chindsR(r),:),hanning(N),[],N,fsamp);
                        Pxy = f13(:,zl(z));
                        if nargin>5
                            Pxy = interp1(F,Pxy,F_scale);
                        else
                            Pxy =  Pxy(F>4);
                        end
                        Pxy = Pxy; %.*tukeywin(length(Pxy),0.25)';
                        xcsd(p,r,z,:) = Pxy;
                    end
                end
            end
        end
        meannpd(chI,chJ,:,:) = (mean(mean(xcsd,1),2));
        clear xcsd
    end
end
F = F_scale;
% % take out symetrical CSD
% diaginds = [2 1; 3 1; 3 2];
% for i = 1:3
%     meancsd(diaginds(i,1),diaginds(i,2),:) = zeros(1,size(meancsd,3));
% end

end