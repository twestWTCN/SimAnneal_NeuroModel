function [P_out,pspT] = convSpikePSP(spT,epsp_t,ampJit,t)
        % Convert spike times to train
        pspT = zeros(size(t));
        absspT = unique(spT);
        spReps = histc(spT,absspT);
        pspT(absspT) = spReps+ (spReps*ampJit.*randn(size(absspT)));
        pspT(pspT<0) = 0;
        % Convolve the train with PSP
        P_out = conv(pspT,epsp_t);
        % Truncate
        P_out = P_out(1:length(t));
        
%         P_out(length(t)+1:end) = [];
