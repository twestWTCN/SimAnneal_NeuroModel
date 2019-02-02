function pnew = resampleParameters_240717(R,pold,prec,MN,initflag)
% computes initial draws from normal distribution (diagonal covariance)
xMin = R.SimAn.pOptBound(1); % set bound on draws (not used in practicality unless prior precision extremely wide)
xMax = R.SimAn.pOptBound(2);
plist = R.SimAn.pOptList;
pnew = pold;
for i = 1:length(plist) % over parameters
    if ~isempty(strfind(plist{i},'src'))
        L = MN;
    else
        L = 1;
    end
    for src = 1:L % over sources
        if isequal(plist{i},'.A') || isequal(plist{i},'.B') % Case for matrices where we loop through
            for Ai = 1:2
                X_old = eval(['pold' plist{i} '{Ai}']);
                if initflag % if initializing then use prior precision
                    X_S = eval(['prec' plist{i} '_s{Ai}']);
                else % else use the newly computed precisions
                    X_S = eval(['prec' plist{i} '{Ai}']);
                end
                Xa = zeros(size(X_old,1),1);
                for Q = 1:numel(X_old)
                    Xa(Q,1) = pearsrnd(X_old(Q),X_S(Q),0,3,1,1); % Kurtosis set to 3 (normal)
                end
                Xa(Xa>xMax & Xa>-28) = xMax; Xa(Xa<xMin & Xa>-28) = xMin;
                X = reshape(Xa,size(X_old));
                
                X(X<-30) = -32;
                X(X_S(:)==0) = X_old(X_S(:)==0);
                eval(['pnew' plist{i} '{Ai} = X;']);
                eval(['pnew' plist{i} '_s{Ai} = X_S;']);
                %                 eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
            X_old = eval(['pold' plist{i}]);
            if initflag % if initializing then use prior precision
                X_S = eval(['prec' plist{i} '_s']);
            else % else use the newly computed precisions
                X_S = eval(['prec' plist{i} '']);
            end
            Xa = zeros(size(X_old,1),1);
            for Q = 1:numel(X_old)
                Xa(Q,1) = pearsrnd(X_old(Q),X_S(Q),0,3,1,1); % Kurtosis set to 3 (normal)
            end
            Xa(Xa>xMax & Xa>-28) = xMax; Xa(Xa<xMin & Xa>-28) = xMin;
            X = reshape(Xa,size(X_old));
            X(X<-30) = -32;
            X(X_S(:)==0) = X_old(X_S(:)==0);
            eval(['pnew' plist{i} ' = X;']);
            eval(['pnew' plist{i} '_s = X_S;']);
        end
    end
end