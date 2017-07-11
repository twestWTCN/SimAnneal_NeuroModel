%function [f] = simAn_master_fx_bgc(x,u,P,M)
clear

%% MASTER FX
load('C:\Users\twest\Documents\Work\PhD\LitvakProject\SimAnneal_NeuroModel\Projects\Rat_CSD\priors\rat_inversion\rat_inversion_params.mat')
% Temp inside
nstates = [8 2 2 2 2 2];
obstate = [3 1 1 1 1 1];
ninputs = [7 1 3 1 3 1];
dt = .001;
tend = 20;
nt = tend/dt;
tvec = linspace(0,tend,nt);
clear x u

for src = 1:6
    x{src} = rand(1,sum(nstates(src)));
    %     exin{src} = rand(1,sum(ninputs(src)));
    
end

%%
fxn{1} = @simAn_fx_mmc; % 8 States
fxn{2} = @simAn_fx_bgc_str; % 2
fxn{3} = @simAn_fx_bgc_gpe; % 2
fxn{4} = @simAn_fx_bgc_stn; % 2
fxn{5} = @simAn_fx_bgc_gpi; % 2
fxn{6} = @simAn_fx_bgc_thal; % 2

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
xstore = zeros(tend,6);

for src = 1:6
    u{src} = rand(1,nstates(src)-1)/32;
    if src == 1 % if source is not MMC
        m(src) = M;
        P.int{1}.C = ones(1,nstates(src));
        P.int{1}.A = repmat(0,1,size(src_input{src},1));
        p(src) = P.int{1};
    else
        m(src) = M;
        m(src).BGC_T = M.BGC_T(src-1);
        m(src).BGC_G = M.BGC_G(src_Gind{src});
        P.int{2}.C = ones(1,nstates(src));
        P.int{2}.A = repmat(0,1,size(src_input{src},1));
        p(src) = P.int{2};
        p(src).T = p(src).T(src-1); % Time constants
        p(src).G = p(src).G(src_Gind{src});     
    end
end

for t = 1:nt
    for src = 1:6
        % Find input signals
        extmp = [];
        for i = 1:size(src_input{src},1)
            extmp(i) = x{src_input{src}(i,1)}(src_input{src}(i,2));
        end
        exin{src} = extmp;

        dx = fxn{src}(x{src},u{src},exin{src},p(src),m(src))';
        x{src} = x{src} + (dx*dt);
        xstore(t,src) = x{src}(obstate(src));
    end
    t
end

