function allPD2 = getPhotodiodeTime2(serverRoot,block,win)
pd = readNPY(fullfile(serverRoot,'photodiode.raw.npy'));
ta1 = readNPY(fullfile(serverRoot,'photodiode.timestamps_Timeline.npy'));
tt = tsToT(ta1, numel(pd)); 
[allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
flipsUp = flipsUp(flipsUp >= win(1) & flipsUp<= win(2));
flipsUp = flipsUp(1:end-1);%%
ntrials = numel(block.events.endTrialValues);
allPD2 = flipsUp(1:ntrials);