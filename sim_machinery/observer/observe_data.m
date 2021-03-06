function [xsims_c R wflag] = observe_data(xstore,m,p,R)
for condsel = 1:numel(R.condnames)
    wflag = 0;
    xsims = xstore{condsel}(R.obs.outstates,:);
    % Delete burnin
    if size(xsims,2) > 5*round(R.obs.brn*(1/R.IntP.dt))
        xsims(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
    else
        wflag = 1;
        xsims_c{condsel} = xsims;
        warning('Simulation is shorter than burn length!!!')
        return
    end
    if any(isnan(xsims(:)))
        xsims_c{condsel} = xsims;
        warning('Simulation contains NaNs!!!')
        return
    end
    tvec_obs = R.IntP.tvec;
    tvec_obs(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
    R.IntP.tvec_obs = tvec_obs;
    
    for i = 1:length(R.obs.gainmeth)
        switch R.obs.gainmeth{i}
            case 'leadfield'
                LF = R.obs.LF.*exp(p.obs.LF);
                LFF = zeros(m.m);
                LFF(eye(size(LFF))~=0) = LF;
                xsims = LFF*xsims;
            case 'difference'
                xsims = [xsims(:,1) diff(xsims,1,2)];
            case 'unitvar'
                for j = 1:size(xsims,1)
                    xsims(j,:) = (xsims(j,:) - mean(xsims(j,:)))./std(xsims(j,:));
                end
            case 'mixing'
                %% REPLACE WITH DISTANCE MATRIX
                mixdeg = R.obs.mixing(1).*exp(p.obs.mixing(1));
                sigmix = repmat(1-mixdeg,m.m,1).*eye(m.m);
                sigmix = sigmix + (repmat(mixdeg/(m.m-1),m.m,1).*~eye(m.m));
                xsims = sigmix*xsims;
            case 'submixing'
                %% REPLACE WITH DISTANCE MATRIX
                mixdeg = R.obs.mixing.*exp(p.obs.mixing);
                cortmix = [mixdeg(1) repmat(mixdeg(2),1,m.m-1)];
                cortmix = [1-(mixdeg(1)*(m.m-1)) repmat(mixdeg(1),1,m.m-1)];
                submix = ~eye(m.m-1,m.m).*repmat(mixdeg(1),m.m-1,m.m);
                submix(logical(eye(size(submix)))) = 1-(mixdeg(2));
                
                mix = [cortmix; circshift(submix,1,2)];
                %             m.m = 6;
                %             dm = eye(m.m) + (~eye(m.m).*repmat(mixdeg(2),m.m,m.m));
                %             dm(2:m.m,1) = repmat(mixdeg(1),1,m.m-1);
                %             dm(1,2:m.m) = repmat(mixdeg(1),m.m-1,1);
                %             dm(logical(eye(size(dm)))) = 0;
                %             dm = dm.*0.001
                xsims = mix*xsims;
            case 'lowpass'
                for i = 1:size(xsims,1)
                    x = xsims(i,:);
                    xsims(i,:) = filtfilt(R.obs.lowpass.fwts,1,x);
                end
            case 'highpass'
                for i = 1:size(xsims,1)
                    x = xsims(i,:);
                    xsims(i,:) = filtfilt(R.obs.highpass.fwts,1,x);
                end
            case 'boring'
%                 figure(100)
%                 plot(xsims'); shg
%                 
                acfcheck = []; acffft = []; montoncheck = [];
                for j = 1:size(xsims,1)
                    swX = slideWindow(xsims(j,:), floor(size(xsims(j,:),2)/3), 0);
                    for swin = 1:size(swX,2)
                        A = swX(:,swin);
                        [acf,lags,bounds] = autocorr(A,1000);
                        acfcheck(j,swin) = any(abs(acf(50:end))>0.90);
                        Env = abs(hilbert(A));
                        [acf,lags,bounds] = autocorr(Env,1000);
                        acfEnvcheck(j,swin) = any(abs(acf(50:end))>0.90);
                        Ar = abs(fft(acf));
                        acffft(j,swin) =  any(Ar(5:end-5)>200);
                        montoncheck(j,swin) = mean(diff(Env(1:0.5/R.IntP.dt:end))>0)>0.70;
                    end
                    
                    %                     Xstab(i) = std(diff(abs(hilbert(xsims(i,:)))))<0.005;
                end
                
                if sum(acffft(:))>4
                    disp('SimFx has boring envelope!')
                    %                     close all
                    %                     error('SimFx is Stable (boring)!')
                    wflag = 1;
                elseif sum(acfcheck(:))>=(numel(acfcheck)/2)
                    disp('SimFx is perfectly periodic!!')
                    wflag = 1;
                elseif sum(acfEnvcheck(:))>4 
                    disp('SimFx has periodic envelope!!')
                    wflag = 1;
                elseif sum(montoncheck(:))>4 %any(acfcheck) ||
                    disp('SimFx seems unstable!')
                    %                     close all
                    %                     error('SimFx is Stable (boring)!')
                    wflag = 1;
                    %                 elseif any(Xstab(:))
                    %                     disp('SimFx seems regular!')
                    %                     wflag = 1;
                end
                if wflag == 0
                    a = 1;
                end
%                 pause(2)
        end
    end
    xsims_c{condsel} = xsims;
end
% if R.obs.norm
%     % Normalise
%     for i = 1:size(xstore,1)
%         xtmp = xstore(i,:);
%         xtmp = (xtmp-mean(xtmp))/std(xtmp);
%         %         xtmp = xtmp.*(m(i).lfpgain*exp(p(i).lfpgain)); % Multiply by gainfield
%         xsims(i,:) = xtmp;
%     end
% end

% Amplify to gain
