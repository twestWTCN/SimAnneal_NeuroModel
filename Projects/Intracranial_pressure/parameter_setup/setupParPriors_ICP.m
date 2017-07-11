function P = setupParPriors_ICP()
% state variables which will change over time dfdt
P.Pic   = 0;        % mmHg ICP
P.Ca    = 0;       % ml/mmHg Compliance
P.vs    = 0;        % mmHg dural venous pressure
P.I     = 0;          % ml/s CSF infusion rate

% constant parameters
P.Ro    = 1;      % mmHg x s/ml CSF outflow resistance
P.Rpv   = 0;       % mmHg x s/ml Proximal venous resistance
P.Rf    = 0;   % mmHg s s/ml CSF fluid formation resistance
P.DCa1  = 0;       % ml/mmHg Sigmoid curve amplitude
P.DCa2  = 0;      % ml/mmHg Sigmoid curve amplitude 
P.Can   = 0;       % ADD HERE
P.Kr    = 0;   % mmHg3 x s /ml;
P.KE    = -1.5;       % 
P.t     = 0;         % seconds
P.Qn    = 0;       % ml/s
P.G     = -1;        % ml x /mmHg x 100% /CBF change
P.Pc    = 1;

% R.SimAn.pOptList = {'.Pic','.Ca','.KE','.Ro','.Can','.G','.Pc'};
