function P = setupParPriors_MVAR_orig(m)
P.n = 0;
P.params(:,:,1) = zeros(m.n);
                      
P.params(:,:,2) = zeros(m.n);
                        
P.noisecov      = zeros(m.n);
