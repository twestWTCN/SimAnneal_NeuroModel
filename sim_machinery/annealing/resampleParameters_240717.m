function pnew = resampleParameters_240717(R,pold,stdev,MN)
xMin = R.SimAn.pOptBound(1);
xMax = R.SimAn.pOptBound(2);
plist = R.SimAn.pOptList;
pnew = pold;
for i = 1:length(plist)
    if ~isempty(strfind(plist{i},'src'))
        L = MN;
    else
        L = 1;
    end
    for src = 1:L
        if strmatch(plist{i},'.A')  % Case for matrices where we loop through
            for Ai = 1:2
                X_old = eval(['pold' plist{i} '{Ai}']);
                X_s = eval(['pold' plist{i} '_s{Ai}']);
                Xa = zeros(size(X_old,1),1);
                for Q = 1:numel(X_old)
                    Xa(Q,1) = pearsrnd(X_old(Q),X_s(Q)*stdev,0,3,1,1); % Kurtosis set to 3 (normal)
                end
                Xa(Xa>xMax & Xa>-30) = xMax; Xa(Xa<xMin & Xa>-30) = xMin;
                X = reshape(Xa,size(X_old));
                X(X_old==-32) = -32;
                eval(['pnew' plist{i} '{Ai} = X;']);
                %                 eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
            X_old = eval(['pold' plist{i}]);
            X_S = eval(['pold' plist{i} '_s']);
            Xa = zeros(size(X_old,1),1);
            for Q = 1:numel(X_old)
                Xa(Q,1) = pearsrnd(X_old(Q),X_S(Q)*stdev,0,3,1,1); % Kurtosis set to 3 (normal)
            end
            Xa(Xa>xMax & Xa>-30) = xMax; Xa(Xa<xMin & Xa>-30) = xMin;
            X = reshape(Xa,size(X_old));
            X(X_old==-32) = -32;
            eval(['pnew' plist{i} ' = X;']);
            eval(['pnew' plist{i} ' = X;']);
        end
    end
end