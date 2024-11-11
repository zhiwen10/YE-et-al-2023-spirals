function [trace2d,traceAmp,tracePhase,rawTrace] = spiralPhaseMap5(U,dV,t,params,freq,rate,mimg)
lowpass = params.lowpass;
gsmooth = params.gsmooth;
% rate = params.rate;
%%
nSV = size(U,3);
Fs = (1/median(diff(t)))*(1/rate);
x = size(U,1); y = size(U,2);
%%
if nargin ==7
    U = U./mimg;
end
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));
%% trace re-construction
meanTrace = Ur*dV;
meanTrace = double(meanTrace);
tsize = size(meanTrace,2);
%%
if rate ~=  1
    tq = 1:rate:tsize;
    meanTrace = interp1(1:tsize,meanTrace',tq);
    meanTrace = meanTrace';
    tsize = numel(tq);
end
%% filter 2-8Hz
meanTrace = meanTrace -mean(meanTrace ,2);
% filter and hilbert transform work on each column
meanTrace = meanTrace';
[f1,f2] = butter(2, 0.5/(Fs/2), 'high');
rawTrace = filtfilt(f1,f2,meanTrace);
%%
[f1,f2] = butter(2, freq/(Fs/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
%%
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
%%
tracePhase = reshape(tracePhase,tsize,x,y);
traceAmp = reshape(traceAmp,tsize,x,y);
trace2d = reshape(meanTrace,tsize,x,y);
rawTrace = reshape(rawTrace,tsize,x,y);
%%
tracePhase = permute(tracePhase,[2,3,1]);
traceAmp = permute(traceAmp,[2,3,1]);
trace2d = permute(trace2d,[2,3,1]);
rawTrace = permute(rawTrace,[2,3,1]);
end