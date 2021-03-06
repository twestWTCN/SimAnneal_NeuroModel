function pnew = resampleParameters_mVar_110717(R,pold,PSS,PPS,stdev,MN)
xMin = R.SimAn.pOptBound(1);
xMax = R.SimAn.pOptBound(2);
plist = R.SimAn.pOptList;
pnew = pold;
for i = 1:length(plist)
    for src = 1:length(MN)
        if strmatch(plist{i},'.A')  % Case for matrices where we loop through
            for Ai = 1:2
                X_old = eval(['pold' plist{i} '{Ai*2}']);
                X_skew = eval(['PSS' plist{i} '{Ai*2}']);
                X_prec = eval(['PPS' plist{i} '{Ai*2}']);
                Xa = zeros(size(X_old,1),1);
                for Q = 1:numel(X_old)
                    Xa(Q,1) = pearsrnd(X_old(Q),X_prec(Q)*stdev,X_skew(Q),3,1,1); % Kurtosis set to 3 (normal)
                end
                Xa(Xa<-20) = -32;
                Xa(Xa>xMax & Xa>-30) = xMax; Xa(Xa<xMin & Xa>-30) = xMin;
                X = reshape(Xa,size(X_old));
                eval(['pnew' plist{i} '{Ai*2} = X']);
                eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
            X_old = eval(['pold' plist{i}]);
            X_skew = eval(['PSS' plist{i}]);
            X_prec = eval(['PPS' plist{i}]);
            Xa = zeros(size(X_old,1),1);
            for Q = 1:numel(X_old)
                Xa(Q,1) = pearsrnd(X_old(Q),X_prec(Q)*stdev,X_skew(Q),3,1,1); % Kurtosis set to 3 (normal)
            end
            Xa(Xa<-20) = -32;
            Xa(Xa>xMax & Xa>-30) = xMax; Xa(Xa<xMin & Xa>-30) = xMin;
            X = reshape(Xa,size(X_old));
            eval(['pnew' plist{i} ' = X;']);
            eval(['pnew' plist{i} ' = X;']);
        end
    end
end