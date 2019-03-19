function epsp_t = makeEPSP(win,tau1,tau2)
epsp_t = (exp(-win/tau2) - exp(-win/tau1));

epsp_t2 = (epsp_t-min(epsp_t))./(max(epsp_t)-min(epsp_t));