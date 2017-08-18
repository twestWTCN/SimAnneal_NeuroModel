function [permMod] = modelProbs(x,m,p,R,d)
if nargin<5
    d = sprintf('%d',[R.d(1:3)]);
end
load([R.rootn 'outputs\' R.out.tag '\parBank_' R.out.tag '_' d '.mat'])
parOptBank = varo;
% figure
% hist(parOptBank(end,:),[-1:.1:1]); xlim([-1 1])
eps = R.analysis.modEvi.eps;
N = R.analysis.modEvi.N;
% parOptBank = parOptBank(:,parOptBank(end,:)>eps);
%% Resample parameters
% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m); % in structure form
pIndMap = spm_vec(pInd); % in flat form
R.SimAn.minRank = ceil(size(pIndMap,1)*1.1);
xf = zeros(size(pIndMap,1),size(parOptBank,2));
for i = 1:size(pIndMap,1)
    x = parOptBank(pIndMap(i),:); % choose row of parameter values
    xf(i,:) = x;
end

disp('Drawing from copula...')
r = copularnd('t',R.Mfit.Rho,R.Mfit.nu,N);
clear x1
for Q = 1:size(xf,1)
    x1(Q,:) = ksdensity(xf(Q,:),r(:,Q),'function','icdf');
end
% setup pars from base
clear base
base = repmat(spm_vec(p),1,N);
for i = 1:N
    base(pIndMap,i) = x1(:,i);
end
spm_unvec(mean(base,2),p)

if isempty(gcp)
    parpool
end

ppm = ParforProgMon('Model Probability Calculation',N);

parfor jj = 1:N
%     ppm.increment();
    pnew = par{jj};
    %% Simulate New Data
        u = innovate_timeseries(R,m);
        u = u./R.IntP.dt;
    xsims = R.IntP.intFx(R,m.x,u,pnew,m);
    % Run Observer function
    if isfield(R.obs,'obsFx')
        xsims = R.obs.obsFx(xsims,m,pnew,R);
    end
    % Run Data Transform
    if isfield(R.obs,'transFx')
        [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,R.obs.SimOrd,R);
    else
        feat_sim = xsims; % else take raw time series
    end
    % Compare Pseudodata with Real
    r2mean  = R.IntP.compFx(R,feat_sim);

%     R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
    r2rep{jj} = r2mean;
    par_rep{jj} = pnew;
    feat_rep{jj} = feat_sim;
    disp(jj); % 
    ppm.increment(); 
end
permMod.r2rep = r2rep;
permMod.par_rep = par_rep;
permMod.feat_rep = feat_rep;
save([R.rootn 'outputs\' R.out.tag '\permMod_' R.out.tag '_' d '.mat'],'permMod')

figure
r2bank = [r2rep{:}];
[h r] = hist(r2bank,25); %D is your data and 140 is number of bins.
h = h/sum(h); % normalize to unit length. Sum of h now will be 1.
bar(h, 'DisplayName', 'Model NRMSE'); 
xD = r(2:2:end);
xL = 2:2:length(r); % list of indices
set(gca,'XTick',xL)
set(gca,'XTickLabel',strread(num2str(xD,2),'%s'))

legend('show');
ylabel('P(D-D*)'); xlabel('D-D*');
hold on
Yval = get(gca,'YLim')

tmp = abs(xD-eps);
[idx idx] = min(tmp); %index of closest value
epsm = xL(xD==xD(idx)); %closest value

plot([epsm epsm],Yval,'B--','linewidth',3)

Pmod = numel(r2bank(r2bank>eps))/N;
annotation(gcf,'textbox',...
    [0.28 0.81 0.19 0.09],...
    'String',{sprintf('eps = %.2f',eps),sprintf('P(m|D) = %.2f',Pmod)},...
    'HorizontalAlignment','right',...
    'FitBoxToText','off',...
    'LineStyle','none');

saveallfiguresFIL_n([R.rootn '\outputs\' R.out.tag '\modelEvidence.jpg'],'-jpg',1,'-r100',1);
