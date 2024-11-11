function [wheel_onset] = getWheelOnsetTime(serverRoot,block)
% get photodiode time
win = [0,5000];
allPD2 = getPhotodiodeTime(serverRoot,win);
ntrial = numel(block.events.endTrialValues);
allPD2 = allPD2(1:ntrial);
%% rotaryEncoder
sigName = 'rotaryEncoder';
tlFile = fullfile(serverRoot, [sigName '.raw.npy']); 
pd = readNPY(tlFile);
tlFile = fullfile(serverRoot, [sigName '.timestamps_Timeline.npy']);
tlTimes = readNPY(tlFile);
tt = tsToT(tlTimes, numel(pd));  
wh = correctCounterDiscont(pd);
%% wheel onset time
trialWin = [0,5];
wheel_onset = nan(size(allPD2));
trialWin3 = trialWin(1):1/200:trialWin(2);
for i = 1:numel(allPD2)
    clear mon mof
    wt = allPD2(i)+trialWin3;
    wval = interp1(tt,wh,wt);
    [mon, mof] = findWheelMoves3(wval, wt, 200);
    if not(isempty(mon))
        wheel_onset(i,1) = mon(1)-allPD2(i);
    end
end