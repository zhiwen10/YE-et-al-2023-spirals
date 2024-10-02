function [Tall,Ta,Tna,thred,traceA,traceNa] = getAmpEpochs(amp,epochLength,traceraw)
%% bin mean amp for each 2 seconds epochs
amp_mean = mean(amp,1);
epochs = floor(size(amp,2)./epochLength);
for i = 1:epochs
    firstFrameAll(i,1) = 1+(i-1)*epochLength;
    lastFrameAll(i,1) = i*epochLength;
    ampEpoch(i,1) = mean(amp_mean(1,firstFrameAll(i,1):lastFrameAll(i,1)));
end
firstFrame = firstFrameAll;
lastFrame = lastFrameAll;
amp = ampEpoch;
Tall = table(firstFrame,lastFrame,amp);
%%
thred = prctile(ampEpoch,80);
%%
clear firstFrame lastFrame amp
amp = ampEpoch(ampEpoch>=thred);
firstFrame = firstFrameAll(ampEpoch>=thred,1);
lastFrame = lastFrameAll(ampEpoch>=thred,1);
alpha_epochN = size(firstFrame,1);
for i = 1:7
    for kk = 1:alpha_epochN
        traceA(i,kk,:) = traceraw(i,firstFrame(kk):lastFrame(kk));
    end
end
Ta = table(firstFrame,lastFrame,amp);
%%    
clear firstFrame lastFrame amp
amp = ampEpoch(ampEpoch<thred);
firstFrame = firstFrameAll(ampEpoch<thred,1);
lastFrame = lastFrameAll(ampEpoch<thred,1);
nonalpha_epochN = size(firstFrame,1);
for i = 1:7
    for kk = 1:nonalpha_epochN
        traceNa(i,kk,:) = traceraw(i,firstFrame(kk):lastFrame(kk));
    end
end
Tna = table(firstFrame,lastFrame,amp);