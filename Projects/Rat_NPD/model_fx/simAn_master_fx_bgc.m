function [xt xobs] = simAn_master_fx_bgc(x,u,p,m,dt,t)
%% MASTER FX
src_input = m.src_input;
for src = 1:6
    x{src} = rand(1,sum(m(src).nstates(src)));
    % Find input signals
    extmp = [];
    for i = 1:size(src_input{src},1)
        extmp(i) = x{src_input{src}(i,1)}(src_input{src}(i,2));
    end
    exin{src} = extmp;
    
    dx = m(src).fxn{src}(x{src},u{src}(t,:),exin{src},p(src),m(src))';
    xt{src} = x{src} + (dx*dt);
    xobs(src) = x{src}(m(src).obstate(src));
end


