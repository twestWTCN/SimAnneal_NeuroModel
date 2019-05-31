function TL = defineBurstTimeLockEpoch(BB,TL,cond)
% EpochTime
Bo = 0;
preBo = [Bo(1)+ floor((TL.periodT(1)/1e3)*BB.fsamp):Bo(1)-1]; %pre burst onset
postBo = [Bo(1): Bo(1) + floor((TL.periodT(2)/1e3)*BB.fsamp) + 1]; % post burst onset
epochdef = [preBo(1):postBo(end)];
TL.epochT = linspace(TL.periodT(1),TL.periodT(2),size(epochdef,2));

% LocalEps
localeps = BB.epsAmpfull;
segInds = BB.segInds{cond};
% Work first with lengths
clear BEpoch REpoch PLVpoch dRPdEpoch meanPLV maxAmp minAmp epsCross usedinds
for i = 1:numel(segInds)
    Bo = segInds{i};
    preBo = [Bo(1)+ floor((TL.periodT(1)/1e3)*BB.fsamp):Bo(1)-1]; %pre burst onset
    postBo = [Bo(1): Bo(1) + floor((TL.periodT(2)/1e3)*BB.fsamp) + 1]; % post burst onset
    epochdef = [preBo(1):postBo(end)];
    
    % Convert from full time to SW time
    if preBo(1)>0 && postBo(end)<size(BB.AEnv{cond},2)
        % Find onset Time aligned to beta onset
        A = BB.AEnv{cond}(:,epochdef);%.*hanning(numel(epochdef))'; % amplitude data
        AH = BB.AEnv{cond}(:,epochdef).*hanning(numel(epochdef))'; % amplitude data
        % Find Crossing times with respect to STN onset
        for L = 1:size(BB.AEnv{cond},1)
            if any(AH(L,:)>localeps(L)) % For finding maximums locally
%                 [dum ec] = find(A(L,:)==max(A(L,:)),1,'first');
                [dum ec] = find(AH(L,:)>localeps(L),1,'first');
                epsCross(L) =  TL.epochT(ec);
%                 epsCross(L) = min(BB.Tvec{cond}(BB.segInds{cond}{L}))
            else
                epsCross(L) = NaN;
            end
        end
        TL.onsetT{cond}(:,i) = epsCross;
        A = (A-min(A,2))./std(A,[],2);
        A = A.^2; %.^2;
        TL.amp{cond}(:,:,i) = A;
        raw = BB.data{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % amplitude data
        TL.raw{cond}(:,:,i) =raw;
        BP = BB.BP{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % amplitude data
        TL.BP{cond}(:,:,i) = BP;

        % STR_M2 RP
%         uwRP = unwrap(BB.Phi{cond}(1,:))-unwrap(BB.Phi{cond}(2,:));
%         RP = wrapToPi(uwRP(epochdef));
%         TL.STR_M2_RP{cond}(i) = circ_mean(RP');
%         TL.STR_M2_dRP{cond}(:,i) = diff(uwRP(epochdef));

        %         PLVbase = nanmedian(BB.PLV{cond});
        %         [dum T(1)] = min(abs(BB.SWTvec{cond}-BB.TSw(epochdef(1))));
        %         T(2) = T(1) + floor(sum(abs(periodT/1000))/diff(BB.TSw(1:2)));
        %         if epochdef(end)<size(BB.AEnv{cond},2) && epochdef(1) > 0 && T(2)<=size(BB.PLV{cond},2)
        %             BEpoch(:,:,i) = 1*zscore(BB.AEnv{cond}(:,epochdef),0,2).*hanning(numel(epochdef))'; % ch x time x burstN
        %             REpoch(:,:,i) = 1*zscore(BB.Tvec{cond}(:,epochdef),0,2).*hanning(numel(epochdef))';
        % %             BEpoch(:,:,i) = 0.2*BB.A{cond}(:,epochdef); %.*hanning(numel(epochdef))'; % ch x time x burstN
        % %             REpoch(:,:,i) = 0.5*BB.rawTime{cond}(:,epochdef).*hanning(numel(epochdef))';
        %             PLVpoch(:,i) = 100*(BB.PLV{cond}(1,T(1):T(2))-PLVbase)/PLVbase ;
        %             dRPdEpoch(:,i) = dRPdt(epochdef)';
        %             meanPLV(i) = mean(PLVpoch(:,i)); %computePPC(squeeze(BB.Phi([1 4],Bo)));
        %             maxAmp(i) = max(BB.AEnv{cond}(4,Bo));
        %             minAmp(i) = min(BB.AEnv{cond}(4,preBo));
        %             Segpoch(i) = segL(segInds(i));
        %         end
    end
end
