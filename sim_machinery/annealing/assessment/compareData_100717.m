function r2mean = compareData_100717(R,sim_dat)
switch R.data.datatype
    %% CSD
    case 'CSD'
        NPDemp  = R.data.feat_emp; % empirical
        NPDsim  = sim_dat; % simulated
        for i = 1:size(NPDemp,1)
            for j = 1:size(NPDemp,2)
                switch R.objfx.feattype
                    case 'complex'
                        if i~=j
                            yfx = (squeeze(imag(NPDsim(i,j,:))));
                            ffx = (squeeze(imag(NPDemp(i,j,:))));
                            r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                            
                            yfx = squeeze(real(NPDsim(i,j,:)));
                            ffx = squeeze(real(NPDemp(i,j,:)));
                            r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                            
                            r2loop(i,j) = mean(r);
                            
                        else
                            yfx = squeeze(real(NPDsim(i,j,:)));
                            ffx = squeeze(real(NPDemp(i,j,:)));
                            r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                            
                            r2loop(i,j) = r(2);
                        end
                    case 'absolute'
                        yfx = squeeze(abs(NPDsim(i,j,:)));
                        ffx = squeeze(abs(NPDemp(i,j,:)));
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
        %% NPD
    case 'NPD'
        NPDemp  = R.data.feat_emp; % empirical
        NPDsim  = sim_dat; % simulated
        for i = 1:size(NPDemp,1)
            for j = 1:size(NPDemp,2)
                switch R.objfx.feattype
                    case 'ForRev'
                        if i==j
                            yfx = (squeeze(NPDsim(i,j,1,:)));
                            ffx = (squeeze(NPDemp(i,j,1,:)));
                            r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                            r2loop(i,j) = r(2);
                        else
                            yfx = (squeeze(NPDsim(i,j,2,:)));
                            ffx = (squeeze(NPDemp(i,j,2,:)));
                            r(1) = goodnessOfFit(yfx,ffx,'NRMSE');
                            
                            yfx = (squeeze(NPDsim(i,j,3,:)));
                            ffx = (squeeze(NPDemp(i,j,3,:)));
                            r(2) = goodnessOfFit(yfx,ffx,'NRMSE');
                            
                            r2loop(i,j) = mean(r);
                        end
                end
            end
        end
        switch R.objfx.specspec;
            case 'auto'
                r2mean = mean(diag(r2loop));
                %         r2mean = sum(diag(r2loop));
                
            case 'cross'
                r2mean = mean(r2loop(triu(r2loop)~=0));
                %         r2mean = sum(r2loop(triu(r2loop)~=0));
        end
        %% TIME
    case 'time' % time courses
        TCemp  = R.data.feat_emp; % empirical
        TCsim  = sim_dat; % simulated
        
        for i = 1:size(TCemp,1)
            try
                yfx = TCsim(i,:);
                ffx = TCemp(i,:);
                if size(yfx,1)<size(yfx,2)
                    yfx = yfx';
                end
                if size(ffx,1)<size(ffx,2)
                    ffx = ffx';
                end
                
                
                r = goodnessOfFit(yfx,ffx,'NRMSE');
                r2loop(i) = r;
            catch
                r2loop(i) = -32;
            end
        end
        r2mean = mean(r2loop);
end