function xstore = stepintegrator_delay_Vec(x,m,p,R,u)

xstore = zeros(m.n,R.IntP.nt);
% Compute X inds (removes need to spm_unvec which is slow)
xinds = zeros(size(m.x,2),2);
for i = 1:size(m.x,2)
    if i == 1
    xinds(i,1) = 1;
    xinds(i,2) = size(m.x{i},2);
    else
        xinds(i,1) = xinds(i-1,2)+1;
        xinds(i,2) = xinds(i,1) + (size(m.x{i},2)-1);
    end
end
m.xinds = xinds;

% Integrate
for t = R.IntP.buffer+1:R.IntP.nt
    [xint] = spm_fx_gen2_delay2(x,u(t,:),p,m,R.IntP.dt,xstore(:,t-R.IntP.buffer:t));
    xstore(:,t) = xint;
    x = spm_unvec(xint,m.x);
    if sum(isnan(xint))>0 || sum(imag(xint))>0
        xstore = NaN(m.n,R.IntP.nt);
        disp('Integration broked :(')
        break
    end
end