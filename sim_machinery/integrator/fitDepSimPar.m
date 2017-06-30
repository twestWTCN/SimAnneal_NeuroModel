function R = fitDepSimPar(R,fitscore)
% RMSE = linspace(0,1,50);
% df = -log(RMSE)+0.5;
% df(RMSE<0.1) = min(df(RMSE<0.1))
% plot(RMSE,df)
% shg
if fitscore<0.1
    R.obs.csd.df = R.obs.csd.df;
else
    R.obs.csd.df = -log(fitscore)+0.25;
end

fsamp = 1/R.IntP.dt;    % sample rate
N = floor(fsamp/R.obs.csd.df);  % segment length
R.IntP.tend = (N*R.obs.csd.reps)/fsamp; % Target time

disp(sprintf('The new simulation length is %.2f seconds',R.IntP.tend))
disp(sprintf('The new simulation df is %.2f Hz',R.obs.csd.df))

if fitscore<0.75 && fitscore>=0.5 
    R.SimAn.rep = 24;
elseif fitscore<0.85 && fitscore>=0.75
    R.SimAn.rep = 16;
elseif fitscore<0.9 && fitscore>=0.85
    R.SimAn.rep = 8;
end
disp(sprintf('Now using %f replicates per temperature',R.SimAn.rep))
