function J =findJacobian(R,x,uc,pc,m)
uc = {zeros(size(uc{1}))};
% computes the Jacobian of a function
n=size(x,1);
xvec = x;
R.IntP.nt = R.IntP.buffer + 6*R.IntP.buffer;
fx = R.IntP.intFx(R,x,uc,pc,m);
fx = fx{1}(:,end);
eps=1.e-8;  % could be made better
xperturb= xvec;
for i=1:n
    xperturb(i,:)=xperturb(i,:)+eps;
    fx_st = R.IntP.intFx(R,xperturb,uc,pc,m);
    fx_st = fx_st{1}(:,end);
    J(:,i)=(fx_st-fx)/eps;
    xperturb(i,:)= xvec(i,:);
end