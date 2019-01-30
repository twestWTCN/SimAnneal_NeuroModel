function [permMod, xsimMod] = modelProbs(x,m,p,R)
if ~isfield(R.analysis.BAA,'flag')
    R.analysis.BAA.flag = 0;
end
% load([R.rootn 'outputs\' R.out.tag '\parBank_' R.out.tag '_' d '.mat'])
parOptBank = R.parOptBank;
% figure
% hist(parOptBank(end,:),[-1:.1:1]); xlim([-1 1])
eps = R.analysis.modEvi.eps;

% parOptBank = parOptBank(:,parOptBank(end,:)>eps);
%% Compute KL Divergence
[KL DKL] = KLDiv(R,p,m,parOptBank);
N = R.analysis.modEvi.N;

if R.analysis.BAA.flag == 0

%% Resample parameters
% Compute indices of optimised parameter
pInd = parOptInds_110817(R,p,m.m,2); % in structure form
pIndMap = spm_vec(pInd); % in flat form
% pIndMap (71 x 1)
% xf (71 x size(parOptBank,2) )
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
    base(pIndMap,i) = x1(:,i); % Insert changing pars into whole base pars
    par{i} = spm_unvec(base(:,i),p);
end
plotDistChange_KS(R.Mfit.Rho,R.Mfit.nu,xf,p,pInd,R,1)

    %% Setup Parfor
    % if isempty(gcp)
    %     parpool
    % end
    a = gcp;
    ppm = ParforProgMon('Model Probability Calculation',N);
    parforArg = a.NumWorkers;
elseif R.analysis.BAA.flag
    base = parOptBank(1:end,:);
    % If doing BAA analysis of model
    parforArg =0;
    ppm = [];
    N = 1;
    
    switch R.analysis.BAA.redmeth
        case 'average'
            % Take the expected parameters from distribution
            par = [];
            par{1} = spm_unvec(mean(base,2),p);
        case 'best'
            par = [];
            par{1} = spm_unvec(base(:,end),p);
        case 'UQ'
            nmrse = base(end-1,:);
            X = base(:,nmrse>prctile(nmrse,75));
            par = [];
            par{1} = spm_unvec(mean(X,2),p);
    end
end

%%Plot Example
figure(5)
pnew = par{1};
[r2mean,pnew,feat_sim,wflag] = computeSimData(R,m,pnew,0,1);
wfstr = ones(1,N);
while wfstr(end)>0
    parfor (jj = 1:N, parforArg)
%         for jj = 1:N
        
        %     ppm.increment();
        pnew = par{jj};
        [r2mean,pnew,feat_sim,xsims,wflag] = computeSimData(R,m,pnew,0);
        
        %     R.plot.outFeatFx({},{feat_sim},R.data.feat_xscale,R,1)
        wfstr(jj) = wflag;
        r2rep{jj} = r2mean;
        par_rep{jj} = pnew;
        feat_rep{jj} = feat_sim;
        disp(jj); %
        if ~R.analysis.BAA.flag
            ppm.increment();
            xsims_rep{jj} = [];
        else
            xsims_rep{jj} = xsims;
        end
    end
    if ~R.analysis.BAA.flag
        wfstr(end) = 0;
    end
end

permMod.r2rep = r2rep;
permMod.par_rep = par_rep;
permMod.feat_rep = feat_rep;
permMod.DKL = DKL;
permMod.KL = KL;
xsimMod = xsims_rep;

% mkdir([R.rootn 'outputs\' R.out.tag '2\'])
% save([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')
% load([R.rootn 'outputs\' R.out.tag '2\permMod_' R.out.tag '_' d '.mat'],'permMod')

% Do plotting
if ~R.analysis.BAA.flag
    figure
    r2bank = [permMod.r2rep{:}];
    [h r] = hist(r2bank,50); %D is your data and 140 is number of bins.
    h = h/sum(h); % normalize to unit length. Sum of h now will be 1.
    bar(h, 'DisplayName', 'Model NRMSE');
    xD = r(2:2:end);
    xL = 2:2:length(r); % list of indices
    set(gca,'XTick',xL)
    set(gca,'XTickLabel',strsplit(num2str(xD,2),' '))
    
    legend('show');
    ylabel('P(D-D*)'); xlabel('D-D*');
    hold on
    Yval = get(gca,'YLim');
    
    tmp = abs(xD-eps);
    [idx idx] = min(tmp); %index of closest value
    epsm = xL(xD==xD(idx)); %closest value
    
    plot([epsm epsm],Yval,'B--','linewidth',3)
    
    Pmod = numel(r2bank(r2bank>eps))/R.analysis.modEvi.N;
    annotation(gcf,'textbox',...
        [0.28 0.81 0.19 0.09],...
        'String',{sprintf('eps = %.2f',eps),sprintf('P(m|D) = %.2f',Pmod)},...
        'HorizontalAlignment','right',...
        'FitBoxToText','off',...
        'LineStyle','none');
    set(gcf,'Position',[680 437 1070 541])
    % saveallfiguresFIL_n([R.rootn 'outputs\' R.out.tag '\modelEvidence.jpg'],'-jpg',1,'-r200',1);
end
