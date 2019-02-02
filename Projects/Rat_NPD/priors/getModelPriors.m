function pQ = getModelPriors(m)
for i = 1:m.m
    switch m.dipfit.model(i).source
        case 'GPE'
            % GPe
            pQ(i).G  = [2]*200;   % synaptic connection strengths
            pQ(i).T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
        case 'STN'
            % GPe
            pQ(i).G  = [2]*200;   % synaptic connection strengths
            pQ(i).T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
        case 'THAL'
            pQ(i).G  = [2]*200;   % synaptic connection strengths
            pQ(i).T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
        case 'STR'
            pQ(i).G  = [2]*200;   % synaptic connection strengths
            pQ(i).T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
        case 'MMC'
            pQ(i).G  = [2 4 2 2 2 2 2 2 2 2 4 2 2 2]*200;         % intrinsic connections
            pQ(i).T  = [3 2 12 18];                               % synaptic time constants [mp sp ii dp]
        case 'GPI'
            pQ(i).G  = [2]*200;   % synaptic connection strengths
            pQ(i).T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
    end
    % Convert to seconds
    pQ(i).T = pQ(i).T./1000;
end
