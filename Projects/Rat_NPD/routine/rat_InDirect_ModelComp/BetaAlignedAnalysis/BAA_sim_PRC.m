function [R] = BAA_sim_PRC(R,MP,simtime)
% Comopute simulations by sweeping across data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
% MP = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);

m = MP.m;
P = MP.p;
% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
fsamp = 1/R.IntP.dt;
conStren = [0.0001 1 1.3];
for i = 1:numel(conStren)
    R.obs.brn = 3; % temporarily!
    % Setup pulses for PRC computation
    delayPulse = (1.*fsamp);
    pulseWid = (0.002.*fsamp);
    pulseAmp = 0.1;
    delayJit = 0.2*fsamp;
    pU = zeros(size(R.IntP.tvec));
    pulseStart = [];
    pulseStart(1) = floor(R.obs.brn.*fsamp);
    cnt = 1;
    while (pulseStart(end)+fsamp) < size(R.IntP.tvec,2)
        pU(pulseStart(cnt):pulseStart(cnt)+pulseWid-1) = ones(pulseWid,1);
        pulseStart(cnt+1) = pulseStart(cnt) + delayPulse + round(delayJit.*randn(1));
        cnt = cnt+1;
    end
    pU = pU.*pulseAmp;
    uc_ip{1} = uc;
    uc_ip{2} =  uc_ip{1};
    uc_ip{2}{1}(:,4) = uc{1}(:,4) + pU';
    
    % Setup the simulations
    Pbase = P;
    Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*conStren(i)); %
    R.obs.brn = 0; % temporarily!
    [r2mean,pnew,feat_sim,dum,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
    [r2mean,pnew,feat_sim,dum,xsim_ip{2}]  = computeSimData(R,m,uc_ip{2},Pbase,0);
    
    XL = [];    
    XL(1,:) = xsim_ip{1}{1}(4,:);
    XL(2,:) = xsim_ip{2}{1}(4,:);
    
    XL = ft_preproc_bandpassfilter(XL,1./R.IntP.dt,[12 28],[],'but');
   
    
    % Phases
    phiL = [];
    phiL(1,:) = angle(hilbert(XL(1,:)));
    phiL(2,:) = angle(hilbert(XL(2,:)));
    
    % phase of impulse
    impPhi = phiL(1,pulseStart);
    % phase difference 3 cycles down the line
    dPhi = wrapToPi(unwrap(phiL(1,:))-unwrap(phiL(2,:)));
    impdPhi = dPhi(pulseStart+(0.03.*fsamp));
    % Amplitudes
    AL(1,:) = abs(hilbert(XL(1,:)));
    AL(2,:) = abs(hilbert(XL(2,:)));
    
    dA = AL(1,:)-AL(2,:);
    impdA = dA(pulseStart);
    
    subplot(2,1,2)
    scatter(impPhi,impdPhi)
    hold on
    xlabel('Impulse Phase \phi');
    ylabel('Change in phase \Delta\phi')
    
    subplot(2,1,1)
    scatter(impPhi,impdA)
    hold on
    xlabel('Impulse Phase \phi');
    ylabel('Change in Amplitude')
end
