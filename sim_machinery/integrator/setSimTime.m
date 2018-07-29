function R = setSimTime(R,N)
R.obs.csd.df = 0.75;
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.dt = .0005;
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
dfact = fsamp/(2*2^(R.obs.SimOrd));
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

disp(sprintf('The Simulation is %.2f seconds long',R.IntP.tend))

disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));
