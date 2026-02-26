function allPD1 = getPhotodiodeTime(session_root,win)
% pd = readNPY(fullfile(session_root,'photodiode.raw.npy'));
% tlTimes = readNPY(fullfile(session_root,'photodiode.timestamps_Timeline.npy'));
sigName = 'photodiode';
load(fullfile(session_root,[sigName '_raw.mat']));                         % load pd
load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));         % load tlTimes
tt = tsToT(tlTimes, numel(pd)); 
%%
[allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
flipsUp = flipsUp(flipsUp >= win(1) & flipsUp<= win(2));
flipsUp = flipsUp(1:end-1);
%%
dff_allPD = diff(allPD);
allPD1 = allPD;
indx = (dff_allPD < 0.6);
indx = [1;indx];
allPD1(logical(indx)) = [];