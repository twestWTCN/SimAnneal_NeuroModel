function r2mean = compareCSD_240517(CSDdata,CSDsim,R)
for i = 1:size(CSDdata,1)
    for j = 1:size(CSDdata,2)
        switch R.objfx.feattype
            case 'complex'
                if i~=j
                    yfx = (squeeze(imag(CSDsim(i,j,:))));
                    ffx = (squeeze(imag(CSDdata(i,j,:))));
                    r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                    
                    yfx = squeeze(real(CSDsim(i,j,:)));
                    ffx = squeeze(real(CSDdata(i,j,:)));
                    r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                    
                    r2loop(i,j) = mean(r);
                    
                else
                    yfx = squeeze(real(CSDsim(i,j,:)));
                    ffx = squeeze(real(CSDdata(i,j,:)));
                    r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                    
                    r2loop(i,j) = r(2);
                end
            case 'absolute'
                yfx = squeeze(abs(CSDsim(i,j,:)));
                ffx = squeeze(abs(CSDdata(i,j,:)));
                r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                r2loop(i,j) = r(1);
        end
    end
end
% r2loop = triu(r2loop);
% r2loop = diag(r2loop);
% r2loop = r2loop(1,1);
% r2loop = reshape(r2loop,1,[]);
% r2loop(isnan(r2loop)) = [];
% r2mean = mean(r2loop,2);
switch R.objfx.specspec;
    case 'auto'
        r2mean = mean(diag(r2loop));
%         r2mean = sum(diag(r2loop));
        
    case 'cross'
          r2mean = mean(r2loop(triu(r2loop)~=0));
%         r2mean = sum(r2loop(triu(r2loop)~=0));
        
end
