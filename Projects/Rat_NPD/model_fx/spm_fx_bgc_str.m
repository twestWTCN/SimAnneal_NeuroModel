function [f,J,Q] = spm_fx_bgc_str(x,u,P,M)
% state equations for a neural mass model of the basal ganglia circuit
% models the circuit between striatum, gpe, stn, gpi, and thalamus as a
% single source (no extrinsic connections)
%
% order           cells     states
% 1 = striatum  - ii        x(1,1:2)
%
% G(1,1) = str -> str (-ve self)

% check if intrinsic connections are free parameters
%--------------------------------------------------------------------------
try, P.G; catch, P.G = 0; end

% get dimensions and configure state variables
%--------------------------------------------------------------------------
% x  = spm_unvec(x,M.x);               % neuronal states
n  = size(x,1);                      % number of sources

% [default] fixed parameters
%--------------------------------------------------------------------------
G  = [2]*200;   % synaptic connection strengths
T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
R  = 2/3;                       % slope of sigmoid activation function
% NB for more pronounced state dependent transfer functions use R  = 3/2;

if isfield(M,'BGC_str_G'); G = M.BGC_str_G; end
if isfield(M,'BGC_str_T'); T = M.BGC_str_T; end


% [specified] fixed parameters
%--------------------------------------------------------------------------
if isfield(M,'pF')
    try, E = M.pF.E; end
    try, G = M.pF.G; end
    try, T = M.pF.T; end
    try, R = M.pF.R; end
end

% input connections
%--------------------------------------------------------------------------
C    = exp(P.C);

% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = R.*exp(P.S);              % gain of activation function
F    = 1./(1 + exp(-R*x + 0));   % firing rate
S    = F - 1/(1 + exp(0));       % deviation from baseline firing (0)

% input
%==========================================================================
% if isfield(M,'u')
% 
%     % endogenous input
%     %----------------------------------------------------------------------
%     U = u(:)*2048;
% 
% else
%     % exogenous input
%     %----------------------------------------------------------------------
%     U = C*u(:)*32;
% end

U = C*u(:)*32;
% time constants and intrinsic connections
%==========================================================================
T     = T/1000;
for i = 1:size(P.T,2)
    T(:,i) = T(:,i).*exp(P.T(:,i));
end


% intrinsic/extrinsic connections to be optimised
%--------------------------------------------------------------------------
j     = 1;
for i = 1:size(P.G,2)
    G(:,j(i)) = G(:,j(i)).*exp(P.G(:,i));
end


% Motion of states: f(x)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - Str: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  U;
u      =  - G(:,1)*S(:,1) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(1,1))./T(1,1);

% Voltage
%==========================================================================
f(:,1) = x(:,2);
f = f';
% f      = spm_vec(f);

% order           cells     states
% 1 = striatum  - ii        x(1,1:2)
if nargout < 2; return, end

% Jacobian
%==========================================================================
if isfield(M,'x'), x = spm_vec(M.x); else,  x = sparse(M.n,1); end
if isfield(M,'u'), u = spm_vec(M.u); else,  u = sparse(M.m,1); end
J  = spm_diff(M.f,x,u,P,M,1);


if nargout < 3; return, end


% delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_dcm_delay(P,M,J);
 

return

