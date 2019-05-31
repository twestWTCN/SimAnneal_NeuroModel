function PRC = computePRC(PRC,XL,pulseStart,bdef,fsamp,cond)
% Computes PRC outputs using difference between perturbed (XL(2xN)) and
% unperturbed (XL(1xN)) time series. A list of pulseStart (1xnPulses)
% indices is given alongside sample rate (fsamp).
    XL = ft_preproc_bandpassfilter(XL,fsamp,bdef,[],'but');
    
    PRC.pulseStart_cS{cond} = pulseStart;
    
    % Phases
    phiL = [];
    phiL(1,:) = angle(hilbert(XL(1,:)));
    phiL(2,:) = angle(hilbert(XL(2,:)));
    PRC.phiL_cS{cond} = phiL;
    
    % phase of impulse
    impPhi = phiL(1,pulseStart);
    PRC.impPhi_cS{cond} = impPhi;
    
    % phase difference 3 cycles down the line
    dPhi = wrapToPi(unwrap(phiL(1,:))-unwrap(phiL(2,:)));
    impdPhi = dPhi(pulseStart+(0.015.*fsamp));
    PRC.impdPhi_cS{cond} = impdPhi;
    
    % Amplitudes
    AL =[];
    AL(1,:) = abs(hilbert(XL(1,:)));
    AL(2,:) = abs(hilbert(XL(2,:)));
    PRC.AL_cS{cond} = AL;
    
    dA = AL(1,:)-AL(2,:);
    impdA = dA(pulseStart);
    PRC.impdA_cS{cond} = impdA;