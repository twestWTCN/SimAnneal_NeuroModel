function pInd = parOptInds_110817(R,p,MN)
plist = R.SimAn.pOptList;
for i = 1:length(plist)
    if ~isempty(strfind(plist{i},'src'))
        L = MN;
    else
        L = 1;
    end
    for src = 1:L
        if strmatch(plist{i},'.A')  % Case for matrices where we loop through
            for Ai = 1:2
                X= eval(['p' plist{i} '{Ai}']);
                x = reshape(X,1,[]);
                xseq = rand(size(x));
                eval(['p' plist{i} '{Ai} = xseq;']);
                pvec = full(spm_vec(p));
                K = strfind(pvec',xseq);
                indK = K:K+(size(x,2)-1);
                indK(x==-32) = [];
                eval(['pInd' plist{i} '{Ai} = indK;']);
%                 eval(['pnew' plist{i} '{(Ai*2)-1} = X']);
            end
            
        else
                X = eval(['p' plist{i}]);
                x = reshape(X,1,[]);
                xseq = rand(size(x));
                eval(['p' plist{i} '= xseq;']);
                pvec = full(spm_vec(p));
                K = strfind(pvec',xseq);
                indK = K:K+(size(xseq,2)-1);
                indK(x==-32) = [];
                eval(['pInd' plist{i} ' = indK;']);
        end
    end
end




% vX = []
% f   = fieldnames(p);
% for i = 1:numel(f)
%     j = 1;
%     for j = 1:numel(p.(f{i}))
%         if iscell(p.(f{i}))
%             X = reshape(p.(f{i}){j},[],1);
%             for k = 0:size(X,1)-1
%             vX = [vX; X repmat(i,size(X)) repmat(j,size(X))];
%             vfield{i,j+k} = ['p.' f{i} '{' num2str(j) '}']
%             end
%         end
%         if isnumeric(p.(f{i}))
%             X = p.(f{i})(j)
%             vX = [vX; X i j];
%             vfield{i,j} = ['p.' f{i} '(' num2str(j) ')']
%         end
%         if isstruct(p.(f{i}))
%         end
%         
%         
%     end
%     
% end
%     X = spm_vec({p.(f{i})});
%     vX = [vX; X repmat(i,size(X))];
% end
% pvec = full(spm_vec(p));
% 
% for i = 1:length(plist)
%     ip = strmatch(plist{i}(2:3),f)
%         for j = 1:numel(p.(f{7}))
%         
% 
% 
% 