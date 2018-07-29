function M = setupParVals_ICP()
% state variables which will change over time dfdt
M.Pic   = 9.5;        % mmHg ICP
M.Ca    = 0.15;       % ml/mmHg Compliance

% input variables
a              = repmat(100,1,2000);
b              = repmat(90,1,2500);
c              = repmat(90,1,10*1e3);
M.Pa           = [a c b];                                                                                   % mmHg Arterial Pressure
M.T            = [1:length(a), length(a)+1e-3:1e-3:length(a)+10, length(a)+11:length(a)+10+length(b)];      % small temporal step size for integration
M.T = linspace(M.T(1),M.T(end),length(M.T));
M.dPadt        = [0 diff(M.Pa)];

    

M.vs    = 6.0;        % mmHg dural venous pressure
M.I     = 0;          % ml/s CSF infusion rate

% constant parameters
M.Ro    = 6*526.3;      % mmHg x s/ml CSF outflow resistance
M.Rpv   = 1.24;       % mmHg x s/ml Proximal venous resistance
M.Rf    = 2.38*1e3;   % mmHg s s/ml CSF fluid formation resistance
M.DCa1  = 0.75;       % ml/mmHg Sigmoid curve amplitude
M.DCa2  = 0.075;      % ml/mmHg Sigmoid curve amplitude 
M.Can   = 0.15;       % ADD HERE
M.Kr    = 4.91*1e4;   % mmHg3 x s /ml;
M.KE    = [ 2.1 0.44 0.27];       % 
M.t     = 20;         % seconds
M.Qn    = 12.5;       % ml/s
M.G     = 1.5;        % ml x /mmHg x 100% /CBF change
M.Pc    = 25;
