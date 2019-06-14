function PRC = computePRC(PRC,XL,pulseStart,bdef,fsamp,cond)
% Computes PRC outputs using difference between perturbed (XL(2xN)) and
% unperturbed (XL(1xN)) time series. A list of pulseStart (1xnPulses)
% indices is given alongside sample rate (fsamp).
    XL_M2 = ft_preproc_bandpassfilter(squeeze(XL(:,1,:)),fsamp,bdef,[],'but');
    XL_STN = ft_preproc_bandpassfilter(squeeze(XL(:,2,:)),fsamp,bdef,[],'but');
    
    PRC.pulseStart_cS{cond} = pulseStart;
    
    % Relative Phases
    phiL = [];
    phiL(1,:) = unwrap(angle(hilbert(XL_STN(1,:)))) - unwrap(angle(hilbert(XL_M2(1,:))));
    phiL(2,:) = unwrap(angle(hilbert(XL_STN(2,:)))) - unwrap(angle(hilbert(XL_M2(2,:))));
    RP = wrapToPi(phiL);
    PRC.phiL_cS{cond} = RP;
    
    % phase of impulse
    impPhi = RP(1,pulseStart);
    PRC.impPhi_cS{cond} = impPhi;
    
    % phase difference 2 cycles down the line
    dPhi = wrapToPi(unwrap(phiL(1,:))-unwrap(phiL(2,:)));
    impdPhi = dPhi(pulseStart+(0.01.*fsamp));
    PRC.impdPhi_cS{cond} = impdPhi;
    
    % Amplitudes
    AL =[];
    AL(1,:) = abs(hilbert(XL_STN(1,:)));
    AL(2,:) = abs(hilbert(XL_STN(2,:)));
    PRC.AL_cS{cond} = AL;
    
    dA = AL(1,:)-AL(2,:);
    impdA = dA(pulseStart+(0.01.*fsamp));
    PRC.impdA_cS{cond} = impdA;