function parsweep = modelBetaParSweep(m,p,parsweep,R)
R.obs.gainmeth = R.obs.gainmeth(1:2);
% R.obs.gainmeth = {'unitvar','leadfield','submixing'};
bwid = 1.5;

u = innovate_timeseries(R,m);
u = u./R.IntP.dt;

N = 12;
Qlist = linspace(-2,2,N);
Rlist = linspace(-2,2,N);
for q = 1:N
    for r = 1:N
        pa = p;
        eval(['pa' parsweep.Q ' = Qlist(q);'])
        eval(['pa' parsweep.R ' = Rlist(r);'])
        psweep{q,r} = pa;
    end
end
betaPowBank = zeros(m.m,N,N);
frqPowBank = cell(N,N);
for q = 1:N
    parfor r = 1:N
        frqlist = 12:2:24;
        pnew = psweep{q,r};
        pnew.A{1}(3,4)
        % Integrate
        xsims = R.IntP.intFx(R,m.x,u,pnew,m);
        % Run Observer function
        if isfield(R.obs,'obsFx')
            xsims = R.obs.obsFx(xsims,m,pnew,R);
        end
        
        % Plot
        % figure
        % plotRawTimeSeries(R,xsims)
        % set(gcf,'Position',[680  281  824  696])
        frqPow = zeros(m.m,length(frqlist));
        for i = 1:m.m
            for j = 1:length(frqlist)
                frqPow(i,j) = bandpower(xsims(i,:),1/R.IntP.dt,[frqlist(j)-bwid frqlist(j)+bwid]);
            end
        end
        betaPow = sum(frqPow(:,frqlist>=14 & frqlist<22),2);
        
        frqPowBank{q,r} = frqPow
        betaPowBank(:,q,r) = betaPow
        disp([q r])
    end
end

parsweep.Rlist = Rlist;
parsweep.Qlist = Qlist;
parsweep.R = '.A{1}(3,4)';
parsweep.Q = '.A{2}(4,3)';
parsweep.frqPowBank = frqPowBank;
parsweep.betaPowBank = betaPowBank;
