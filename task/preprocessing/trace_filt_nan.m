function [traceFilt2,tracePhase2] = trace_filt_nan(trace_mean)
% trace_mean shape: x * y * t
trace_mean2 = reshape(trace_mean,[size(trace_mean,1)*size(trace_mean,2),size(trace_mean,3)]);
%%
Fs = 35;
freq = [2,8];
%% take care of nan pixels after registration before filting 
indx1 = not(isnan(trace_mean2(:,1)));
trace_mean4 = trace_mean2(indx1,:);
meanTrace = trace_mean4 -mean(trace_mean4,2);
meanTrace = meanTrace';  
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
traceFilt = filtfilt(f1,f2,meanTrace);
traceHilbert =hilbert(traceFilt);
tracePhase = angle(traceHilbert);
tracePhase2 = nan(size(trace_mean2));
traceFilt2 = nan(size(trace_mean2));
tracePhase2(indx1,:) = tracePhase';
traceFilt2(indx1,:) = traceFilt';
traceFilt2 = reshape(traceFilt2,[size(trace_mean,1),size(trace_mean,2),size(trace_mean,3)]);
tracePhase2 = reshape(tracePhase2,[size(trace_mean,1),size(trace_mean,2),size(trace_mean,3)]);
