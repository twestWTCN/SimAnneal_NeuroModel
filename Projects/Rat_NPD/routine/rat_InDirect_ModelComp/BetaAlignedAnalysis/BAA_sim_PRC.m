function [R] = BAA_sim_PRC(R,MP,simtime)
% Comopute simulations by sweeping across data
% [R,m,permMod,xsimMod{1}] = getSimModelData_v2(R,modID,simtime);
% MP = permMod{1}.par_rep{1};
R = setSimTime(R,simtime);
cmap = brewermap(18,'Spectral');
PRC.condcmap = cmap([1 4 8 16 4 18],:);
PRC.condname = {'Fitted','1% STN->GPe','150% STN->GPe'}; %Fitted','1% M2->STN','150% M2->STN',
m = MP.m;
P = MP.p;
% connection list
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
fsamp = 1/R.IntP.dt;
conStren = [0.3 1 1.3];
for cond = 1:numel(conStren)
    R.obs.brn = 3; % temporarily!
    % Setup pulses for PRC computation
    delayPulse = (1.*fsamp);
    pulseWid = (0.002.*fsamp);
    pulseAmp = 0.1;
    delayJit = 0.2*fsamp;
    pU = zeros(size(R.IntP.tvec));
    pulseStart = [];
    pulseStart(1) = floor(R.obs.brn.*fsamp);
    cnt = 1;
    while (pulseStart(end)+fsamp) < size(R.IntP.tvec,2)
        pU(pulseStart(cnt):pulseStart(cnt)+pulseWid-1) = ones(pulseWid,1);
        pulseStart(cnt+1) = pulseStart(cnt) + delayPulse + round(delayJit.*randn(1));
        cnt = cnt+1;
    end
    pU = pU.*pulseAmp;
    uc_ip{1} = uc;
    uc_ip{2} =  uc_ip{1};
    uc_ip{2}{1}(:,1) = uc{1}(:,1) + pU'; % Give it a cortical pulse
    
    % Setup the simulations
    Pbase = P;
    Pbase.A{1}(3,4) = log(exp(Pbase.A{1}(3,4))*conStren(cond)); %
    R.obs.brn = 0; % temporarily!
    [~,~,~,~,xsim_ip{1}] = computeSimData(R,m,uc_ip{1},Pbase,0);
    [~,~,~,~,xsim_ip{2}]  = computeSimData(R,m,uc_ip{2},Pbase,0);
    
    XL = [];
    XL(1,:,:) = xsim_ip{1}{1}([1 4],:);
    XL(2,:,:) = xsim_ip{2}{1}([1 4],:);
    
    PRC = computePRC(PRC,XL,pulseStart,[12 28],fsamp,cond);
end
save([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\PRCtmp'],'PRC')
load([R.rootn '\routine\' R.out.tag '\BetaBurstAnalysis\Data\PRCtmp'],'PRC')

plotSimulatedPRC(PRC,[2 1 3])
set(gcf,'Position',[711         604        1081         374])