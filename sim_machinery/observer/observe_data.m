function [xsims R] = observe_data(xstore,m,p,R)

xsims = xstore([7 9 11 13 15 17],:);
% Delete burnin
xsims(:,1:round(R.obs.brn*(1/R.IntP.dt))) = [];
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
        case 'unitvar'
            for i = 1:size(xsims,1)
                xsims(i,:) = (xsims(i,:) - mean(xsims(i,:)))./std(xsims(i,:));
            end
        case 'mixing'
                mixdeg = R.obs.mixing.*exp(p.obs.mixing);
                sigmix = repmat(1-mixdeg,m.m,1).*eye(m.m);                
                sigmix = sigmix + (repmat(mixdeg/(m.m-1),m.m,1).*~eye(m.m));
                xsims = sigmix*xsims;
    end
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
