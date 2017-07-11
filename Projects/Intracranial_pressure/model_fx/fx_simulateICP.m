function [xint,T] = fx_simulateICP(x,m,p)
% Log Normal Priors
a = fieldnames(p);
P = m;
for i = 1:numel(a)
    y = eval(['m.' a{i}]); ys = eval(['p.' a{i}]);
    Y = y.*exp(ys);
    eval(['P.' a{i} '= Y;'])
end

for j = 1%:length(P.KE)
    P.Ke = P.KE(j);
    try
        [T, f]  = ode15s(@(t,F) odefx(t,F,m.Z,P),P.T,x);
        VolArt  = f(:,2)'.*(P.Pa - f(:,1)');
        Ra      = [];
        xint = f';
    catch
        f = NaN;
        xint = NaN(size(f));
    end
    %     for t = 1:length(T)
%         [~,r] = odefx(T(t),f(t,:),Z,P);
%         Ra    = [Ra r];
%     end
%     CBF     = (P.Pa - P.Pc)./Ra;
    
    
end

function [dfdt,Ra] = odefx(t,F,Z,P)
% initiate equations of motion

% Z(1)    = P.Pa;
% Z(2)    = P.dPadt;
% Z(3)    = P.vs;
% Z(4)    = P.I;

% F(1)    = P.Pic;
% F(2)    = P.Ca;

% dfdt(1) = dCadt;
% dfdt(2) = dPicdt;

z1      = interp1(P.T,Z{1},t);
z2      = interp1(P.T,Z{2},t);


Va      = F(2)*(z1 - F(1));
Ra      = P.Kr * (P.Can^2)/(Va^2);
Pc      = ((z1*P.Rpv) + (F(1)*Ra))/(P.Rpv +Ra);
q       = (z1 - Pc)/Ra;
x       = (q - P.Qn)/P.Qn;
if x<0
   DCa  = P.DCa1;   
else
   DCa  = P.DCa2;
end
Ksig    = DCa/4;

SigGx   = ((P.Can + DCa/2) + (P.Can - DCa/2)*exp(P.G*x/Ksig)) / (1+ exp(P.G*x/Ksig));
%SigGx   = Sig(P,x,Ksig,DCa);


dCadt   = (1/P.t)*(-F(2) + SigGx);

t1 = (P.Ke*F(1))/(1 + (F(2)*P.Ke*F(1)));
t2 = F(2)*z2;
t3 = dCadt*(z1 - F(1));
t4 = (Pc - F(1))/P.Rf;
t5 = (F(1)- Z{3})/P.Ro;

dPicdt  = t1*(t2+t3+t4-t5+Z{4});
dfdt    = [dPicdt;dCadt];

% 
% function s = Sig(P,x,Ksig,DCa)
% 
% s = ((P.Can + DCa/2) + (P.Can - DCa/2)*exp(P.G*x/Ksig)) / (1+ exp(P.G*x/Ksig));
