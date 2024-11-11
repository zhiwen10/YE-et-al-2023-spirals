%% flow
function [N,edges,phase_mu,phase_var,flow_mu,flow_var] = get_flow_metric1(traceAmp_raw,tracePhase1_raw,tracePhase1_pred,vxy_raw,vxy_predict)
amp_all = squeeze(sum(traceAmp_raw,[2,3],'omitnan'))./sum(sum(not(isnan(squeeze(traceAmp_raw(1,:,:)))))); % mean amp
edges = [0:0.00125:0.025];
% [N,edges] = histcounts(amp_all,10);
[N,edges] = histcounts(amp_all,edges);
for i = 1:numel(edges)-1
    clear a alpha beta flow_diff
    a = find(amp_all>=edges(i) & amp_all<edges(i+1));
    indx1{i} = a;
    %%
    if not(isempty(a))
        
        phase_alpha = tracePhase1_raw(indx1{i},:,:);
        phase_alpha(isnan(phase_alpha(:))) = [];
        phase_beta = tracePhase1_pred(indx1{i},:,:);
        phase_beta(isnan(phase_beta(:))) = [];
        %%
        phase_diff = angdiff(phase_alpha(:),phase_beta(:));
        [phase_var(i)] = circ_var(phase_diff);
        [phase_mu(i)] = circ_mean(phase_diff);
            %%
        flow_alpha = vxy_raw(indx1{i},:,:);
        flow_alpha(isnan(flow_alpha(:))) = [];
        flow_alpha = angle(flow_alpha);
        flow_beta = vxy_predict(indx1{i},:,:);
        flow_beta(isnan(flow_beta(:))) = [];
        flow_beta = angle(flow_beta);
        %%
        flow_diff = angdiff(flow_alpha(:),flow_beta(:));
        [flow_var(i)] = circ_var(flow_diff);
        [flow_mu(i)] = circ_mean(flow_diff);
    else
        phase_mu(i) = nan;
        phase_var(i) = nan;
        flow_mu(i) = nan;
        flow_var(i) = nan;
    end       
end
end