function xstore = stepintegrator_delay_Vec(x,m,p,R,u)
% temp variable
tstep = 1000;
% Initialise the state space
xstore = zeros(m.n,R.IntP.nt);

%% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
    xinds(i,1) = 1;
    xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;


%% Define indexes of sources and sinks
efferent(4,:) = [3 3 7 7];               % sources of MMC connections
% afferent(4,:) = [2 4 8 0];               % targets of MMC connections
afferent(4,:) = [8 2 8 4];                  % forward deep/middle; back deep/superficial

efferent(5,:) = [1 1 1 1];               % sources of STR connections
afferent(5,:) = [2 2 2 2];               % targets of STR connections

efferent(6,:) = [1 1 1 1];               % sources of GPE connections
afferent(6,:) = [2 2 2 2];               % targets of GPE connections

efferent(7,:) = [1 1 1 1];               % sources of STN connections
afferent(7,:) = [2 2 2 2];               % targets of STN connections

efferent(8,:) = [1 1 1 1];               % sources of GPI connections
afferent(8,:) = [2 2 2 2];               % targets of GPI connections

efferent(9,:) = [1 1 1 1];               % sources of THAL connections
afferent(9,:) = [2 2 2 2];               % targets of THAL connections

% get the neural mass models {'ERP','CMC'}
n     = numel(x);
model = m.dipfit.model;
for i = 1:n
    if  strcmp(model(i).source,'ERP')
        nmm(i) = 1;
    elseif strcmp(model(i).source,'CMC')
        nmm(i) = 2;
    elseif strcmp(model(i).source,'BGC')
        nmm(i) = 3; 
     elseif strcmp(model(i).source,'MMC')
        nmm(i) = 4;
    elseif  strcmp(model(i).source,'STR')
        nmm(i) = 5;
    elseif  strcmp(model(i).source,'GPE')
        nmm(i) = 6;
    elseif  strcmp(model(i).source,'STN')
        nmm(i) = 7;    
    elseif  strcmp(model(i).source,'GPI')
        nmm(i) = 8;    
    elseif  strcmp(model(i).source,'THAL')
        nmm(i) = 9;
    end
end

%%

% Compute value of delays from lognormal mean
D = zeros(m.m);
D(p.D>-30) = 4/1000; % set all delay priors to 4ms.

D(4,1) = 10/1000;   % M1 to STR (Gerfen and Wilson 1996)
D(2,1) = 2.5/1000;  % M1 to STN (Nakanishi et al. 1987)
D(3,2) = 2/1000;    % STN to GPe (Kita and Kitai 1991)
D(2,3) = 4/1000;    % GPe to STN (Fujimoto and Kita 1993)
D(5,4) = 5/1000;    % STR to GP (Kita and Kitai 1991)

D(6,5) = 4/1000;    % GPe to Thal
D(1,6) = 6/1000;    % Thal to M1


D = ceil(D.*exp(p.D).*(1/R.IntP.dt)); % As expectation of priors and convert units to steps
D(D>R.IntP.buffer) = R.IntP.buffer-1; % Ensure not bigger than nback - 2nd check in script!!
D(D<((1e-3)/R.IntP.dt)&D>0) = floor((1e-3)/R.IntP.dt); % Minimum 1ms   

Ds = zeros(size(D));
% Now find indices of inputs
% Currently no seperation between inh and excitatory
for i = 1:length(nmm) % target
    for j = 1:length(D(i,:)) % source
        if D(i,j)>0
        Ds(i,j) = afferent(nmm(j),1); % input sources
        Ds(i,j) = (m.xinds(j,1)-1)+Ds(i,j);
        end
    end
end
% Integrate
for t = R.IntP.buffer+1:R.IntP.nt
    [xint] = spm_fx_gen2_delay2(x,u(t,:),p,m,R.IntP.dt,xstore(:,t-R.IntP.buffer:t));
    xstore(:,t) = xint;
    x = spm_unvec(xint,m.x);
    if sum(isnan(xint))>0 || sum(imag(xint))>0
        xstore = NaN(m.n,R.IntP.nt);
        disp('Integration broked :(')
        break
    end
end