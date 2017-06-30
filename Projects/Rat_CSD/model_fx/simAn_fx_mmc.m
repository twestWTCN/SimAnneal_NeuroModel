function [f,J,Q] = simAn_fx_mmc(x,u,exin,P,M)
% state equations for a neural mass model of motor cortex
% Bhatt et al. 2016 Neuroimage
%
% FORMAT [f,J,D] = spm_fx_mmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_mmc(x,u,P,M)
% FORMAT [f]     = spm_fx_mmc(x,u,P,M)
% x      - state vector
%   x(:,1) - voltage     (middle pyramidal cells)
%   x(:,2) - conductance (middle pyramdidal cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%   x(:,7) - voltage     (deep pyramidal cells)
%   x(:,8) - conductance (deep pyramidal cells)
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward, backward, lateral) extrinsic rates 
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% T  = synaptic time constants
% S  = slope of sigmoid activation function
%
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% [default] fixed parameters
%--------------------------------------------------------------------------
G  = [2 4 2 2 2 2 2 2 2 2 4 2 2 2]*200;         % intrinsic connections
T  = [3 2 12 18];                               % synaptic time constants [mp sp ii dp]

G = M.MMC_G;
T = M.MMC_T;

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
% for i = 1:length(A)
%     L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
%     A{i} = A{i}./(1 + 4*L);
% end
 
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = 2/3;                      % slope of sigmoid activation function
B    = 0;                        % bias or background (sigmoid)
R    = R.*exp(P.S);              % gain of activation function
F    = 1./(1 + exp(-R*x + B));   % firing rate
S    = F - 1/(1 + exp(B));       % deviation from baseline firing

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
FEx    = 1./(1 + exp(-R*exin + B));   % firing rate
SEx    = FEx - 1/(1 + exp(B));       % deviation from baseline firing
% input
%==========================================================================
C = exp(P.C);
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*1536;
    
else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*32;
    
end

E    = [1]*200;     % extrinsic (forward and backward)
A = exp(P.A).*E;      % forward  connections (sp -> mp)
% time constants and intrinsic connections
%==========================================================================
T    = T/1000;
G    = G;

% extrinsic connections
%--------------------------------------------------------------------------
% forward  (i)   2  sp -> mp (+ve)
% forward  (ii)  1  sp -> dp (+ve)
% backward (i)   2  dp -> sp (-ve)
% backward (ii)  1  dp -> ii (-ve)
%--------------------------------------------------------------------------
% free parameters on time constants and intrinsic connections
%--------------------------------------------------------------------------
% G(:,1)  mp -> mp (-ve self)  4
% G(:,2)  mp -> sp (+ve rec )  4
% G(:,3)  ii -> mp (-ve rec )  4
% G(:,4)  ii -> ii (-ve self)  4
% G(:,5)  mp -> ii (+ve rec )  4
% G(:,6)  dp -> ii (+ve rec )  2
% G(:,7)  sp -> sp (-ve self)  4
% G(:,8)  sp -> mp (+ve rec )  4
% G(:,9)  ii -> dp (-ve rec )  2
% G(:,10) dp -> dp (-ve self)  1
% G(:,11) sp -> dp (+ve rec)  2
% G(:,12) ii -> sp (-ve rec)  4
% G(:,13) sp -> ii (+ve rec)  4
% G(:,14) dp -> sp (+ve rec)  2
%--------------------------------------------------------------------------
% Neuronal states (deviations from baseline firing)
%--------------------------------------------------------------------------
%   S(:,1) - voltage     (middle pyramidal cells)
%   S(:,2) - conductance (middle pyramidal cells)
%   S(:,3) - voltage     (superficial pyramidal cells)
%   S(:,4) - conductance (superficial pyramidal cells)
%   S(:,5) - current     (inhibitory interneurons)
%   S(:,6) - conductance (inhibitory interneurons)
%   S(:,7) - voltage     (deep pyramidal cells)
%   S(:,8) - conductance (deep pyramidal cells)
%--------------------------------------------------------------------------
    T = T.*exp(P.T);
    G = G.*exp(P.G);

% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
 
% Middle layer (middle pyramidal cells): Hidden causes
%--------------------------------------------------------------------------
u      =   U(1)+ (A*SEx(1)); %A{1}*S(:,3)+ U;
u      = - G(:,1).*S(:,1) - G(:,3).*S(:,5) + G(:,8).*S(:,3) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(:,1))./T(:,1);
 
% Supra-granular layer (superficial pyramidal cells): Hidden causes - error
%--------------------------------------------------------------------------
u      =   A*S(:,3);% + U;
u      =  - G(:,7).*S(:,3) + G(:,2).*S(:,1) - G(:,12).*S(:,5) + G(:,14).*S(:,7) + u;
f(:,4) =  (u - 2*x(:,4) - x(:,3)./T(:,2))./T(:,2);
 
% Supra-granular layer (inhibitory interneurons): Hidden states - error
%--------------------------------------------------------------------------
u      =  - A*S(:,7);
u      =  - G(:,4).*S(:,5) + G(:,5).*S(:,1) + G(:,6).*S(:,7) + G(:,13).*S(:,3) + u;
f(:,6) =  (u - 2*x(:,6) - x(:,5)./T(:,3))./T(:,3);
 
% Infra-granular layer (deep pyramidal cells): Hidden states
%--------------------------------------------------------------------------
u      =   A*S(:,7);% + U;
u      = - G(:,10).*S(:,7) - G(:,9).*S(:,5) + G(:,11).*S(:,3) + u;
f(:,8) =  (u - 2*x(:,8) - x(:,7)./T(:,4))./T(:,4);
 
% Voltage
%==========================================================================
f(:,1) = x(:,2);
f(:,3) = x(:,4);
f(:,5) = x(:,6);
f(:,7) = x(:,8);
f      = f';
 


