function [x,p,m] = setup_new_priors(P,M)


%% SPECIFY MODEL PARAMETERS
% each fx file takes (x,u,P,M) as an argument that is:
% x = state
% u = exogenous noise
% P = parameters
% M = model specification and priors

% order           cells     states
% 1 = striatum  - ii        x(1,1:2)
% 2 = gpe       - ii        x(1,3:4)
% 3 = stn       - pyr       x(1,5:6)
% 4 = gpi       - ii        x(1,7:8)
% 5 = thalamus  - pyr       x(1,9:10)
%
% G(1,1) = str -> str (-ve self)
% G(1,2) = str -> gpe (-ve ext)
% G(1,3) = gpe -> gpe (-ve self)
% G(1,4) = gpe -> stn (-ve ext)
% G(1,5) = stn -> gpe (+ve ext)
% G(1,6) = str -> gpi (-ve ext)
% G(1,7) = stn -> gpi (+ve ext)
% G(1,8) = gpi -> gpi (-ve self)
% G(1,9) = gpi -> tha (-ve ext)

% Indices
M.nstates = [8 2 2 2 2 2]; % number of states
M.obstate = [3 1 1 1 1 1]; % indice of state giving rise to observed
M.dobs = [1 4 3 2]; % indices of data corresponding to sources
M.lfpgain = 1;
M.fxn{1} = @simAn_fx_mmc; % 8 States
M.fxn{2} = @simAn_fx_bgc_str; % 2
M.fxn{3} = @simAn_fx_bgc_gpe; % 2
M.fxn{4} = @simAn_fx_bgc_stn; % 2
M.fxn{5} = @simAn_fx_bgc_gpi; % 2
M.fxn{6} = @simAn_fx_bgc_thal; % 2

src_Gind{2} = [1]; % STR (-ve self)
src_Gind{3} = [2 3 5]; % GPE (-ve ext),(-ve self),(+ve ext)
src_Gind{4} = [4]; % STN (-ve ext)
src_Gind{5} = [6 7 8]; % GPI (-ve ext),(+ve ext),(-ve self)
src_Gind{6} = [9]; % THAL (-ve ext)

% source inputs (source and state)
src_input{1} = [6 1]; % MMC - (Thalamus,Voltage)
src_input{2} = [2 1; 1 7]; % STR - Self,(MMC,Deep)
src_input{3} = [3 1; 2 1; 4 1]; % GPe - Self,(Thalamus,Voltage),(STN,Voltage)
src_input{4} = [3 1; 1 7]; % STN - (GPe,Voltage),(MMC,Deep)
src_input{5} = [5 1; 4 1; 2 1]; % GPI - Self,(STN,Voltage),(STR,Voltage)
src_input{6} = [5 1]; % THAL - (STN,Voltage);

M.src_input = src_input;
for src = 1:6
    if src == 1 % if source is not MMC
        m(src) = M;
        P.int{1}.lfpgain = 0;
        P.int{1}.C = repmat(0,1,m(src).nstates(src));
        P.int{1}.A = repmat(0,1,size(src_input{src},1));
        p(src) = P.int{1};
    else
        m(src) = M;
        m(src).BGC_T = M.BGC_T(src-1);
        m(src).BGC_G = M.BGC_G(src_Gind{src});
        P.int{2}.lfpgain = 0;
        P.int{2}.C = repmat(0,1,m(src).nstates(src));
        P.int{2}.A = repmat(0,1,size(src_input{src},1));
        p(src) = P.int{2};
        p(src).T = p(src).T(src-1); % Time constants
        p(src).G = p(src).G(src_Gind{src});
    end
    x{src} = rand(1,sum(m(src).nstates(src))); % initial conditions
end

