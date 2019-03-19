function P_out = convSpikePSP(spT,epsp_t,t)
        % Convert spike times to train
        pspT = zeros(size(t));
        pspT(spT) = 1+ 0.05.*randn(size(spT));
        % Convolve the train with PSP
        P_out = conv(pspT,epsp_t);
        % Truncate
        P_out = P_out(1:length(t));
%         P_out(length(t)+1:end) = [];
