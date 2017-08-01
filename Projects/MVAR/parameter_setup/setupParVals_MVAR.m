function M = setupParVals_MVAR()
M.n = 3;
M.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];
                      
M.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];
                        
M.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];
M.fn =  fieldnames(M);