function PRC = computePRC(PRC,XL,pulseStart,bdef,fsamp,cond)
% Computes PRC outputs using difference between perturbed (XL(i=2xN)) and
% unperturbed (XL(i=1xN)) time series. A list of pulseStart (1xnPulses)
% indices is given alongside sample rate (fsamp).
    XL_M2 = ft_preproc_bandpassfilter(squeeze(XL(:,1,:)),fsamp,bdef,[],'but');
    XL_STN = ft_preproc_bandpassfilter(squeeze(XL(:,2,:)),fsamp,bdef,[],'but');
    
    PRC.pulseStart_cS{cond}(B_sel) = pulseStart;
    
    % Relative Phases
    phiL_STN = [];
    phiL_STN(1,:) = unwrap(angle(hilbert(XL_STN(1,:))));
    phiL_STN(2,:) = unwrap(angle(hilbert(XL_STN(2,:))));
    
    RPL = [];
    RPL(1,:) = unwrap(angle(hilbert(XL_STN(1,:)))) - unwrap(angle(hilbert(XL_M2(1,:))));
    RPL(2,:) = unwrap(angle(hilbert(XL_STN(2,:)))) - unwrap(angle(hilbert(XL_M2(2,:))));
    RP = wrapToPi(RPL);
    PRC.phiL_cS{cond} = RP;
    
    % phase of impulse
    impPhi = RP(1,pulseStart);
    PRC.impPhi_cS{cond} = impPhi;
    
    % phase difference 2 cycles down the line
    dPhi = wrapToPi(unwrap(phiL_STN(1,:))-unwrap(phiL_STN(2,:)));
    impdPhi = dPhi(pulseStart+(0.015.*fsamp));
    PRC.impdPhi_cS{cond} = impdPhi;
    
    % Amplitudes
    AL =[];
    AL(1,:) = abs(hilbert(XL_STN(1,:)));
    AL(2,:) = abs(hilbert(XL_STN(2,:)));
    PRC.AL_cS{cond} = AL;
    
    bLev =  100.*(AL(1,:)-median( AL(1,:)))./median( AL(1,:));
    PRC.impBetaLev{cond} = bLev(pulseStart);
    
    dA = AL(1,:)-AL(2,:);
    impdA = dA(pulseStart+(0.015.*fsamp));
    PRC.impdA_cS{cond} = impdA;