function [freq1,psdx,psdx_mean] = fft_spectrum(trace)
% trace is ntrace x timestamps
alpha_trace_mean = bsxfun(@minus, trace, mean(trace, 2)); % mean-subtract each SVD
N = size(alpha_trace_mean,2);
xdft = fft(alpha_trace_mean');
Fs = 35;
xdft1 = xdft(1:N/2+1,:);
nf = size(xdft1,1);
psdx = (1/(Fs*N)) * abs(xdft1).^2;
psdx(2:end-1,:) = 2*psdx(2:end-1,:);
%%
freq1 = 0:Fs/N:Fs/2;
psdx_mean = mean(psdx,2);