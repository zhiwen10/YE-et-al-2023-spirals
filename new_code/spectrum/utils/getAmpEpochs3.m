function [Tall,Ta,Tna,thred1,thred2,traceA,traceNa] = getAmpEpochs3(traceAmp,epochLength,rawTrace)
%% bin mean amp for each 2 seconds epochs
amp_mean = mean(traceAmp(3:7,:),1);                                       % only use SSp traces
epochs = floor(size(traceAmp,2)./epochLength);
for i = 1:epochs
    firstFrameAll(i,1) = 1+(i-1)*epochLength;
    lastFrameAll(i,1) = i*epochLength;
    ampEpoch(i,1) = mean(amp_mean(1,firstFrameAll(i,1):lastFrameAll(i,1)));
end
firstFrame = firstFrameAll;
lastFrame = lastFrameAll;
traceAmp = ampEpoch;
Tall = table(firstFrame,lastFrame,traceAmp);
%%
thred1 = prctile(ampEpoch,80);
thred2 = prctile(ampEpoch,20);
%%
clear firstFrame lastFrame traceAmp
traceAmp = ampEpoch(ampEpoch>=thred1);
firstFrame = firstFrameAll(ampEpoch>=thred1,1);
lastFrame = lastFrameAll(ampEpoch>=thred1,1);
alpha_epochN = size(firstFrame,1);
for i = 1:8
    for kk = 1:alpha_epochN
        traceA(i,kk,:) = rawTrace(i,firstFrame(kk):lastFrame(kk));
    end
end
Ta = table(firstFrame,lastFrame,traceAmp);
%%    
clear firstFrame lastFrame traceAmp
traceAmp = ampEpoch(ampEpoch<thred2);
firstFrame = firstFrameAll(ampEpoch<thred2,1);
lastFrame = lastFrameAll(ampEpoch<thred2,1);
nonalpha_epochN = size(firstFrame,1);
for i = 1:8
    for kk = 1:nonalpha_epochN
        traceNa(i,kk,:) = rawTrace(i,firstFrame(kk):lastFrame(kk));
    end
end
Tna = table(firstFrame,lastFrame,traceAmp);