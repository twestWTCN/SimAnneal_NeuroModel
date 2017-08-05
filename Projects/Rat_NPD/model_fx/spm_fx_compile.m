function [xstore] = spm_fx_compile(R,x,u,p,m)
% Compiles NMM functions with delays, seperates intrinsic and extrinsic
% dynamics then summates
xinds = m.xinds;
% model-specific parameters
%==========================================================================

% model or node-specific state equations of motions
%--------------------------------------------------------------------------
fx{1} = @spm_fx_erp;                                    % ERP model
fx{2} = @spm_fx_cmc;                                    % CMC model
fx{3} = @spm_fx_bgc;                                    % basal ganglia circuit
fx{4} = @spm_fx_mmc;                                    % motor micro circuit

fx{5} = @spm_fx_bgc_str;
fx{6} = @spm_fx_bgc_gpe;
fx{7} = @spm_fx_bgc_stn;
fx{8} = @spm_fx_bgc_gpi;
fx{9} = @spm_fx_bgc_thal;
% indices of extrinsically coupled hidden states
%--------------------------------------------------------------------------
efferent(1,:) = [9 9 9 9];               % sources of ERP connections
afferent(1,:) = [4 8 5 8];               % targets of ERP connections

efferent(2,:) = [3 3 7 7];               % sources of CMC connections
afferent(2,:) = [2 8 4 6];               % targets of CMC connections

efferent(3,:) = [9 9 9 9];               % sources of BGC connections (thalamus)
afferent(3,:) = [2 6 2 6];               % targets of BGC connections (striatum & STN)

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



% scaling of afferent extrinsic connectivity (Hz)
%--------------------------------------------------------------------------
E(1,:) = [1 0 1 0]*200;                    % ERP connections
E(2,:) = [1 .3571 1 .625]*100000;          % CMC connections (to ctx) with T = [2 2 16 28] gives [200 100 200 100] = regular DCM
E(3,:) = [1.8 1.2 1.8 1.2]*100000;         % BGC connections (to bgc) with T_str=8 and T_stn=4 gives A = 144 and 48
E(4,:) = [.1111 .6667 1 .1111]*100000;             % MMC connections (to mmc) with T_mp=3 and T_sp=2 gives A = 270 and 180; with T_dp=18 gives A=200

%% to calculate E divide the target value for A by the value of the time constant (in seconds, i.e. 0.018)
E(5,:) = [.5 .5 -.5 -.5]*100000;             % STR connections
E(6,:) = [.5 .5 -.5 -.5]*100000;             % GPE connections
E(7,:) = [ 1  1 -.1  -1]*100000;             % STN connections
E(8,:) = [.5 .5 -.5 -.5]*100000;               % GPI connections
E(9,:) = [.5 .5 -.5 -.5]*100000;               % THAL connections

% E(5,:) = [.5 .5 .5 .5]*100000;             % STR connections
% E(6,:) = [.5 .5 .5 .5]*100000;             % GPE connections
% E(7,:) = [ 1  1 .1  1]*100000;             % STN connections
% E(8,:) = [.5 .5 .5 .5]*100000;               % GPI connections
% E(9,:) = [.5 .5 .5 .5]*100000;               % THAL connections


% get the neural mass models {'ERP','CMC'}
%--------------------------------------------------------------------------
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

%% Pre-integration extrinsic connection parameters

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
D(D<((1e-3)/R.IntP.dt)&D>0) = floor((1e-3)/R.IntP.dt); % Minimum 1ms

if (R.IntP.buffer-max(max(D)))<=0
    R.IntP.buffer = max(max(D)) + 2;
    disp(['Delay is bigger than buffer, increasing buffer to: ' num2str(R.IntP.buffer)])
end
Ds = zeros(size(D));Dt = zeros(size(D));
% Now find indices of inputs
% Currently no seperation between inh and excitatory
for i = 1:length(nmm) % target
    for j = 1:length(D(i,:)) % source
        if D(i,j)>0
            Ds(i,j) = afferent(nmm(j),1); % input sources
            Ds(i,j) = (m.xinds(j,1)-1)+Ds(i,j);
            Dt(i,j) = efferent(nmm(i),1); % input sources
            Dt(i,j) = (m.xinds(i,1)-1)+Dt(i,j);
        end
    end
end

% synaptic activation function
%--------------------------------------------------------------------------
Rz     = 2/3;                      % gain of sigmoid activation function
B     = 0;                        % bias or background (sigmoid)
Rz     = Rz.*exp(p.S);
S     = @(x,Rz,B)1./(1 + exp(-Rz*x(:) + B)) - 1/(1 + exp(B));
% dSdx  = @(x,R,B)(R*exp(B - R*x(:)))./(exp(B - R*x(:)) + 1).^2;
% for i = 1:n
%     Sx{i} = S(x{i},R,B);
% end

% Extrinsic connections
%--------------------------------------------------------------------------
alist = [1 2; 3 4];
for i = 1:numel(p.A)
    A{alist(i,1)} = exp(p.A{i});
    A{alist(i,2)} = exp(p.A{i});
end

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
% TOL   = exp(2);
% for i = 1:numel(A)
%     L    = (A{i} > TOL) & (A{i}' > TOL);
%     A{i} = A{i}./(1 + 4*L);
% end

% and scale of extrinsic connectivity (Hz)
%--------------------------------------------------------------------------
for i = 1:n
    for k = 1:numel(p.A)
        A{k}(i,:) = E(nmm(i),k)*A{k}(i,:);
    end
end

%% TIME INTEGRATION STARTS HERE ===========================================
f = zeros(xinds(end),1); dt = R.IntP.dt;
xstore= full(repmat(spm_vec(x),1,R.IntP.buffer)); xint = ones(m.n,1);
TOL = exp(3);
for tstep = R.IntP.buffer:R.IntP.nt
    % assemble flow
    %==========================================================================
    N     = m;
    for i = 1:n % targets
        % intrinsic flow
        %----------------------------------------------------------------------
        N.x  = m.x{i};
        ui   = u(tstep,i);
        Q    = p.int{i};
        Q.C  = p.C(i,:);
        xi = xstore(m.xinds(i,1):m.xinds(i,2),tstep)';
        f(m.xinds(i,1):m.xinds(i,2)) = fx{nmm(i)}(xi,ui,Q,N);
        % extrinsic flow
        %----------------------------------------------------------------------
        for j = 1:n % sources
            for k = 1:numel(p.A) % connection type
                if abs(A{k}(i,j)) > TOL
                    %                 ik       = afferent(nmm(i),k);
                    %                 jk       = efferent(nmm(j),k);
                    %                 xd = spm_unvec(xback(:,end-D(i,j)),M.x);
                        xD = xstore(Ds(i,j),tstep-D(i,j));
                    
                    f(Dt(i,j)) = f(Dt(i,j)) + A{k}(i,j)*S(xD,Rz,B);
                end
            end
        end
    end
    xint = xint + (f.*dt);
    xstore = [xstore xint];
    % if tstep >R.IntP.buffer*10
    %     pause
    % end
    % disp(tstep/R.IntP.nt)
    % xint= spm_unvec(x,M.x);
end
