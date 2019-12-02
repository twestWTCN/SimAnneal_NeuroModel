function plotTrackPhase(stngpephi)

for cond = 1:2
    CTX_phi = squeeze(stngpephi{cond}(1,:,:));
    GPe_phi = squeeze(stngpephi{cond}(2,:,:));
    STN_phi = squeeze(stngpephi{cond}(3,:,:));
    
    
    
    for seg = 1:size(CTX_phi,2)
        t_list = find(diff(STN_phi(:,seg))<-3);
        RP_CTX_STN_seg = wrapToPi(CTX_phi(:,seg)-STN_phi(:,seg));
        RP_GPe_STN_seg = wrapToPi(GPe_phi(:,seg)-STN_phi(:,seg));
        for i = 1:numel(t_list)
            rpCTXSTN(i,seg) = RP_CTX_STN_seg(t_list(i));
            rpGPeSTN(i,seg) = RP_GPe_STN_seg(t_list(i));
        end
    end
    
    rpCTXSTN(rpCTXSTN==0) = nan;
    plot(nanmean(rpCTXSTN,2))
    hold on
    rpCTXSTN(rpGPeSTN==0) = nan;
    plot(nanmean(rpGPeSTN,2))
end