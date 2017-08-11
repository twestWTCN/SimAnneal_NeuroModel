function pInd = parOptInds_110817(R,p,MN)
plist = R.SimAn.pOptList;
pvec = full(spm_vec(p));
for i = 1:length(plist)
    for src = 1:MN
        if strmatch(plist{i},'.A')  % Case for matrices where we loop through
            for Ai = 1:2
                X= eval(['p' plist{i} '{Ai}']);
                X = reshape(X,1,[]);
                K = strfind(pvec',X);
                indK = K:K+(size(X,2)-1);
                indK(X==-32) = [];
                eval(['pInd' plist{i} '{Ai} = indK;']);
%                 eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
            X = eval(['p' plist{i}]);
                X = reshape(X,1,[]);
                K = strfind(pvec',X);
                indK = K:K+(size(X,2)-1);
                indK(X==-32) = [];
            eval(['pInd' plist{i} ' = indK;']);
        end
    end
end
