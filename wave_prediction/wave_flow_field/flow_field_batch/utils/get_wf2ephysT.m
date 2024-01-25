function [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT(ops,t)
%% transform widefield time to ephys time
syncTL = loadAlign(ops.serverRoot, 'tl',ops.imec);
syncProbe = loadAlign(ops.serverRoot, ops.probeName,ops.imec);
if numel(syncProbe)>numel(syncTL)
    syncProbe = syncProbe(1:length(syncTL));
end
%%
WF2ephysT1 = interp1(syncTL, syncProbe, t); 
% WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));