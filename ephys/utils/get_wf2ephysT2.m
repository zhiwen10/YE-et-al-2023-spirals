function [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t)
%% transform widefield time to ephys time
syncTL = readNPY(fullfile(ops.session_root, 'tl_sync.npy'));
syncProbe = readNPY(fullfile(ops.session_root, [ops.probeName '_sync.npy']));
if numel(syncProbe)>numel(syncTL)
    syncProbe = syncProbe(1:length(syncTL));
end
%%
WF2ephysT1 = interp1(syncTL, syncProbe, t); 