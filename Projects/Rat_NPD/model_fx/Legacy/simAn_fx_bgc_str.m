function [f,J,Q] = simAn_fx_bgc_str(x,u,exin,P,M)
% state equations for a neural mass model of the basal ganglia circuit
% models the circuit between striatum, gpe, stn, gpi, and thalamus as a
% single source (no extrinsic connections)
%
% order           cells     states
% 1 = striatum  - ii        x(1,1:2)
%
% G(1,1) = str -> str (-ve self)



% [default] fixed parameters
%--------------------------------------------------------------------------
G  = [2]*200;   % synaptic connection strengths
T  = [8];               % synaptic time constants [str,gpe,stn,gpi,tha];
R  = 2/3;                       % slope of sigmoid activation function
% NB for more pronounced state dependent transfer functions use R  = 3/2;



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

FEx    = 1./(1 + exp(-R*exin));   % firing rate
SEx    = FEx - 1/(1 + exp(0));       % deviation from baseline firing

E    = [1]*200;     % extrinsic (forward and backward)
A = exp(P.A).*E(1);      % forward  connections (sp -> mp)

% input
%==========================================================================
if isfield(M,'u')

    % endogenous input
    %----------------------------------------------------------------------
    U = u(:)*128;

else
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:)*32;
end


% time constants and intrinsic connections
%==========================================================================
    T = T.*exp(P.T);
    G = G.*exp(P.G);

% Motion of states: f(x)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - Str: ii

% inhibitory interneurons: Hidden states - error
%--------------------------------------------------------------------------
u      =  U(1) + (A(1)*SEx(1));
u      =  - G(:,1)*S(:,1) + u;
f(:,2) =  (u - 2*x(:,2) - x(:,1)./T(1,1))./T(1,1);

% Voltage
%==========================================================================
f(:,1) = x(:,2);
f      = f';
% order           cells     states
% 1 = striatum  - ii        x(1,1:2)