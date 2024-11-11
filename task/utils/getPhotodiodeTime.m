function allPD1 = getPhotodiodeTime(serverRoot,win)
pd = readNPY(fullfile(serverRoot,'photodiode.raw.npy'));
ta1 = readNPY(fullfile(serverRoot,'photodiode.timestamps_Timeline.npy'));
tt = tsToT(ta1, numel(pd)); 
%%
[allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
flipsUp = flipsUp(flipsUp >= win(1) & flipsUp<= win(2));
flipsUp = flipsUp(1:end-1);
%%
dff_allPD = diff(allPD);
allPD1 = allPD;
indx = (dff_allPD<0.6);
indx = [1;indx];
allPD1(logical(indx)) = [];