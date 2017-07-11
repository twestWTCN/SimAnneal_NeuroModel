function P = setupPars_ICP()
% state variables which will change over time dfdt
P.Pic   = 9.5;        % mmHg ICP
P.Ca    = 0.15;       % ml/mmHg Compliance

% input variables
a              = repmat(100,1,2000);
b              = repmat(90,1,2500);
c              = repmat(90,1,10*1e3);
P.Pa           = [a c b];                                                                                   % mmHg Arterial Pressure
P.T            = [1:length(a), length(a)+1e-3:1e-3:length(a)+10, length(a)+11:length(a)+10+length(b)];      % small temporal step size for integration

P.dPadt        = [0 diff(P.Pa)];

    

P.vs    = 6.0;        % mmHg dural venous pressure
P.I     = 0;          % ml/s CSF infusion rate

% constant parameters
P.Ro    = 6*526.3;      % mmHg x s/ml CSF outflow resistance
P.Rpv   = 1.24;       % mmHg x s/ml Proximal venous resistance
P.Rf    = 2.38*1e3;   % mmHg s s/ml CSF fluid formation resistance
P.DCa1  = 0.75;       % ml/mmHg Sigmoid curve amplitude
P.DCa2  = 0.075;      % ml/mmHg Sigmoid curve amplitude 
P.Can   = 0.15;       % ADD HERE
P.Kr    = 4.91*1e4;   % mmHg3 x s /ml;
P.KE    = [ 2.1 0.44 0.27];       % 
P.t     = 20;         % seconds
P.Qn    = 12.5;       % ml/s
P.G     = 1.5;        % ml x /mmHg x 100% /CBF change
P.Pc    = 25;
