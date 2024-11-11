function [trace_mean3,traceFilt3,tracePhase3] = taskTrace_upsample(trace_mean_correct,freq,rate)
%%
meanTrace = reshape(trace_mean_correct,size(trace_mean_correct,1)*size(trace_mean_correct,2),...
size(trace_mean_correct,3));
%%
tsize = size(meanTrace,2);
x = size(trace_mean_correct,1); y = size(trace_mean_correct,2);
Fs = 35/rate;
tq = 1:rate:tsize;
meanTrace = interp1(1:tsize,meanTrace',tq);
meanTrace = meanTrace';
tsize = numel(tq);
%% take care of nan pixels after registration before filting 
indx1 = not(isnan(meanTrace(:,1)));
trace_mean4 = meanTrace(indx1,:);
meanTrace2 = trace_mean4 -mean(trace_mean4,2);
meanTrace2 = meanTrace2';  
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
traceFilt = filtfilt(f1,f2,meanTrace2);
traceHilbert =hilbert(traceFilt);
tracePhase = angle(traceHilbert);
tracePhase2 = nan(size(meanTrace));
traceFilt2 = nan(size(meanTrace));
tracePhase2(indx1,:) = tracePhase';
traceFilt2(indx1,:) = traceFilt';
%%
meanTrace = meanTrace-mean(meanTrace,2);
trace_mean3 = reshape(meanTrace,size(trace_mean_correct,1),size(trace_mean_correct,2),size(meanTrace,2));
traceFilt3 = reshape(traceFilt2,size(trace_mean_correct,1),size(trace_mean_correct,2),size(traceFilt2,2));
tracePhase3 = reshape(tracePhase2,size(trace_mean_correct,1),size(trace_mean_correct,2),size(tracePhase2,2));
end