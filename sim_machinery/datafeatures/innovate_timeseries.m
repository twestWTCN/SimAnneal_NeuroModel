function u  = innovate_timeseries(R,m)

fs = 1/R.IntP.dt;
nyq = fs/2;
fw = R.frqzfull./nyq;

if size(fw,2)>size(fw,1)
    fw = fw';
end

if rem(size(fw,1),2)==0
    fw(end) = []; % Make uneven length
end

switch R.IntP.Utype
    case 'DCM_Str_Innov'
        DCM = m.uset.p;
        [Gu,Gs,Gn,f] = spm_csd_mtf_gu(DCM.Ep,fw);
        
        for i = 1:size(Gu,2)
            F =  [fw; 1]; % Freqencies
            M = [Gu(:,i); 0]; % Magnitudes
            b = firls(480,F,M); % Design Filter
            [h1,w] = freqz(b,1); % Filter freq response
            
            x = randn(1,R.IntP.nt+(2*fs)); % White noise
            x = (x-mean(x))/std(x); % normed
            xr = filter(b,1,x); % Apply filter
            xr = xr(fs*1:end-fs*1); % Remove transients
            xr = (xr-mean(xr))./std(xr); % norm
            
            u(:,i) = xr.*m.uset.p.scale; % save filtered noise
        end
    case 'white_covar'
        u = (sqrtm(m.uset.p.covar)*randn(m.m,R.IntP.nt)).*m.uset.p.scale;
        u = u';
end