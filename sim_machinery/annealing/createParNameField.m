function [pNameStruc] = createParNameField(R,p,MN,set)
% This function will find the indices of the parameters to be optimized.
% Parameters will only be added if there expected values are non neglible
% (i.e. > -32).
if nargin<4
    set =2;
end
plist = R.SimAn.pOptList;
for i = 1:length(plist)
    if ~isempty(strfind(plist{i},'src'))
        L = MN;
    else
        L = 1;
    end
    for src = 1:L
        if isequal(plist{i},'.A') || isequal(plist{i},'.B') || isequal(plist{i},'.D') % Case for matrices where we loop through
            for Ai = 1:size(eval(['p' plist{i}]),2)
                X = eval(['p' plist{i} '{Ai}']);
                Atype = {'(+)','(-)'}
                % Now loop through elements of matrix
                nameMat = [];
                for tar = 1:size(X,1)
                    for src = 1:size(X,2)
                        nameMat(tar,src) = [src*100 + tar*10 + Ai];
%                         nameMat{tar,src} = [R.chsim_name{src} ' to ' R.chsim_name{tar} ' ' Atype{Ai}];
                    end
                end
                eval(['pNameStruc' plist{i} '{Ai} = nameMat;']);               
            end
            
        else
            X = eval(['p' plist{i}]);
            x = reshape(X,1,[]);
            
            S = eval(['p' plist{i} '_s']);
            indS = reshape(S,1,[]);
            
            xseq = rand(size(x));
            eval(['p' plist{i} '= xseq;']);
            pvec = full(spm_vec(p));
            indK = [];
            for pin = 1:size(xseq,2)
                indK(pin) = strfind(pvec',xseq(pin));
            end
            
            xseq = rand(size(x));
            eval(['p' plist{i} '_s = xseq;']);
            pvec = full(spm_vec(p));
            indS = [];
            for pin = 1:size(xseq,2)
                indS(pin) = strfind(pvec',xseq(pin));
            end
            
            parSelInd = find(indS==0 | x<-30);
            if set == 1
                indK(parSelInd) = NaN;
                x(parSelInd) = NaN;
                indS(parSelInd) = NaN;
            elseif set == 2
                indK(parSelInd) = [];
                x(parSelInd) = [];
                indS(parSelInd) = [];
            end
            eval(['pInd' plist{i} ' = indK;']);
            eval(['pMu' plist{i} ' = indK;']);
            eval(['pSig' plist{i} ' = indS;']);
            
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