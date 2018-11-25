delta = 1; %eps = eps+1;
while eps_p-eps > 0.001
    pip = 100;
    A = [1];
    while size(A,2)<(R.SimAn.minRank*delta)
        pip = pip-0.001;
        if pip<0
            break
        end
        eps = prctile(parBank(end,:),pip);
        A = parBank(:,parBank(end,:)>eps);
    end
    parOptBank = A;
    clear A
    delta = delta-0.005;
    if delta<0.25
        break
    end
end
