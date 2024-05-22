function [trace2d,traceAmp,tracePhase] = spiralPhaseMap_fftn_freq(U1,dV1,t,params,freq,rate)
lowpass = params.lowpass;
%%
nSV = size(U1,3);
Fs = (1/median(diff(t)))*(1/rate);
x = size(U1,1); y = size(U1,2);
%%
Ur = reshape(U1, size(U1,1)*size(U1,2), size(U1,3));
meanTrace = Ur*dV1;
meanTrace = double(meanTrace);
tsize = size(meanTrace,2);
%% phase scramble
meanTrace = reshape(meanTrace,x,y,tsize); % data
%%
data_fft = fftn(meanTrace); % 3D fourier transform ??????
phase = angle(data_fft);
indx = randperm(numel(phase(:)));
phase1 = phase(indx);
%%
phase1 = reshape(phase1,size(phase));
mag = abs(data_fft) ;
data_fft_new = mag.*exp(1i*phase1);
meanTrace1 = real(ifftn(data_fft_new)) ; % 3D inverse fourier tranform
meanTrace1 = reshape(meanTrace1,x*y, tsize);
%%
if rate ~=  1
    tq = 1:rate:tsize;
    meanTrace1 = interp1(1:tsize,meanTrace1',tq);
    meanTrace1 = meanTrace1';
    tsize = numel(tq);
end
%%
% meanTrace2 = reshape(meanTrace1,size(phase));
% figure;
% for i =1:10
%     subplot(2,10,i);
%     imagesc(squeeze(meanTrace(:,:,1000+i)));
%     axis image; axis off;
%     subplot(2,10,i+10);
%     imagesc(squeeze(meanTrace2(:,:,1000+i)));
%     axis image; axis off;
% end
%% filter 2-8Hz
meanTrace1 = meanTrace1 -mean(meanTrace1 ,2);
% filter and hilbert transform work on each column
meanTrace1 = meanTrace1';
if lowpass
    [f1,f2] = butter(2, 0.2/(Fs/2), 'low');
else
    [f1,f2] = butter(2, freq/(Fs/2), 'bandpass');
end
%%
meanTrace1 = filtfilt(f1,f2,meanTrace1);
%%
traceHilbert =hilbert(meanTrace1);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
%%
tracePhase = reshape(tracePhase,tsize,x,y);
traceAmp = reshape(traceAmp,tsize,x,y);
trace2d = reshape(meanTrace1,tsize,x,y);
%%
tracePhase = permute(tracePhase,[2,3,1]);
traceAmp = permute(traceAmp,[2,3,1]);
trace2d = permute(trace2d,[2,3,1]);
end