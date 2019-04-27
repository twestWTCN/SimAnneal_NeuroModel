clear; close all
load('ksdense_test.mat')        
[Mfit,cflag] = postEstCopula(parOptBank,Mfit,pIndMap,pOrg);
        par = postDrawCopula(R,Mfit,pOrg,pIndMap,pSigMap,rep);
