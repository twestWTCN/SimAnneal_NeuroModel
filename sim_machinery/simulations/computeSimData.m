function [r2mean,pnew,feat_sim,xsims,wflag] = computeSimData(R,m,pnew,simtime,plotop)
if nargin<5
    plotop = 0;
end
if simtime ~= 0
    R = setSimTime(R,simtime);
end
%% Simulate New Data
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
[xsims,tvec,wflag] = R.IntP.intFx(R,m.x,u,pnew,m);
if wflag == 0
    % Run Observer function
    if isfield(R.obs,'obsFx')
        [xsims R] = R.obs.obsFx(xsims,m,pnew,R);
    end
    % Run Data Transform
    if isfield(R.obs,'transFx')
        [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
    else
        feat_sim = xsims; % else take raw time series
    end
    % Compare Pseudodata with Real
    r2mean  = R.IntP.compFx(R,feat_sim);
    disp('Simulation Success!')
    
    if plotop == 1
        xsims{1} = xsims{1} + linspace(m.m,0,m.m)';
        plot(R.IntP.tvec_obs(1:end),xsims{1}(:,2:end))
        legend(R.chsim_name); xlabel('Time (s)'); ylabel('Amplitude')
        set(gcf,'Position',[705         678        1210         420]);
        xlim([4 5])
    end
    
else
    r2mean = -inf;
    feat_sim = NaN;
    disp('Simulation error!')
end