%function [] = modelProbs(x,m,p,varo,R)
%% Resample parameters
Rho = R_out.Mfit.Rho;
nu = R_out.Mfit.nu;
N = 1000;

eps = 0.25;
parOptBank = parBank(:,parBank(end,:)>eps);

r = copularnd('t',Rho,nu,N);
pvec = spm_vec(p);
pvec(149) = 0.1;
ilist = find(pvec>-32 & pvec~= 0);

clear x1
for Q = 1:size(ilist,1)
    x1(Q,:) = ksdensity(parOptBank(ilist(Q),:),r(:,Q),'function','icdf');
end
clear base
base = spm_vec(p);
for jj = 1:N
    base(ilist) = x1(:,i);
    par{jj} = spm_unvec(base,p);
end
if isempty(gcp)
    parpool
end
% ppm = ParforProgMon('Model Probability Calculation',N);

parfor jj = 1:N
%     ppm.increment();
    pnew = par{jj};
    %% Simulate New Data
    % Integrate in time master fx function
    %         xsims = eval([R.IntP.intFx R.IntP.intFxArgs]);
    xsims = R.IntP.intFx(x,m,pnew,R);
    
    if isfield(R.obs,'obsFx') % Run Observer function
        xsims = R.obs.obsFx(xsims,m,pnew,R);
    end
    if isfield(R.obs,'transFx') % Run Data Transform
        %% Construct CSD and compare to data
        %             fx = R.obs.transFx;
        %             [~,feat_sim] = fx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,10,R);
        [~,feat_sim] = R.obs.transFx(xsims,R.chloc_name,R.chloc_name,1/R.IntP.dt,7,R);
    else
        feat_sim = xsims; % else take raw time series
    end
    % Now using NRMSE
    r2mean  = R.IntP.compFx(R,feat_sim);
    R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
    r2rep{jj} = r2mean;
    par_rep{jj} = pnew;
    disp(jj); % 
end

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
