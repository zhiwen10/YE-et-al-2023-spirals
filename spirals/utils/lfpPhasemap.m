function [traceEphys,tracePhase,traceAmp] = lfpPhasemap(ops, epochT,Fs)
%% load downsampled LFP data
lfp_file = ops.fproc1;
fsize  = get_file_size(lfp_file)/384/2; % size in bytes of raw binary
fidW1  = fopen(lfp_file, 'r'); % open for writing processed data
lfp = fread(fidW1, [384 fsize],'int16'); % write this batch to binary file
fclose(fidW1);
%% 
t2 = 0:1/ops.fs*300:(ops.NT*ops.Nbatch-1)/ops.fs;
%% find index of epochT in lfp
[min1,indx1] = min(abs(t2-epochT(1)));
[min2,indx2] = min(abs(t2-epochT(2)));
%% padding 100 sample before and after epoch for interp1 at tTest1 and tTest2
EphysIndx = indx1:indx2;
EphysIndx2 = indx1-100:indx2+100;
meanTrace = lfp(:,EphysIndx2);
%% interp1 at exact time of tTest1 and tTest2
% Fs = 100;
rate = 1/Fs;
tq1 = epochT(1):rate:epochT(2);
meanTrace1 = interp1(t2(EphysIndx2),meanTrace',tq1);
meanTrace1 = meanTrace1';
%% filter 2-8Hz
traceEphys1 = meanTrace1-repmat(mean(meanTrace1,2),1,size(meanTrace1,2));
% filter and hilbert transform work on each column
traceEphys1 = traceEphys1';
[f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
traceEphys = filtfilt(f1,f2,traceEphys1);
traceHilbert =hilbert(traceEphys);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
end  