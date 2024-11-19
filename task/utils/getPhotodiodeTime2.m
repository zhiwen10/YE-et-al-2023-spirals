function allPD2 = getPhotodiodeTime2(session_root,block,win)
% pd = readNPY(fullfile(serverRoot,'photodiode.raw.npy'));
% tlTimes = readNPY(fullfile(serverRoot,'photodiode.timestamps_Timeline.npy'));
sigName = 'photodiode';
load(fullfile(session_root,[sigName '_raw.mat']));                         % load pd
load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));         % load tlTimes
tt = tsToT(tlTimes, numel(pd));

[allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
flipsUp = flipsUp(flipsUp >= win(1) & flipsUp<= win(2));
flipsUp = flipsUp(1:end-1);%%
ntrials = numel(block.events.endTrialValues);
allPD2 = flipsUp(1:ntrials);