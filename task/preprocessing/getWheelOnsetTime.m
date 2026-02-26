function [wheel_onset] = getWheelOnsetTime(session_root,block)
% get photodiode time
win = [0,5000];
allPD2 = getPhotodiodeTime(session_root,win);
ntrial = numel(block.events.endTrialValues);
allPD2 = allPD2(1:ntrial);
%% rotaryEncoder
sigName = 'rotaryEncoder';
% tlFile = fullfile(session_root, [sigName '.raw.npy']); 
% pd = readNPY(tlFile);
% tlFile = fullfile(session_root, [sigName '.timestamps_Timeline.npy']);
% tlTimes = readNPY(tlFile);
load(fullfile(session_root,[sigName '_raw.mat']));                         % load pd
load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));         % load tlTimes
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