githubDir = 'C:\Users\Steinmetz lab\Documents\git';
% Add paths
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'ephysWidefieldPhasemap')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysWidefieldPhasemap')
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\pinwheelDetection')
%% session info
mn = 'ZYE_0016';
td = '2021-02-25';
en = 1;
imec = 0;
probeName = 'p1';

% mn = 'ZYE_0059';
% td = '2022-03-15';
% en = 1;
% imec = 0;
% probeName = 'p0';

% mn = 'ZYE_0060';
% td = '2022-03-24';
% en = 1;
% imec = 0;
% probeName = 'p0';

% mn = 'ZYE_0062';
% td = '2022-03-30';
% en = 5;
% imec = 0;
% probeName = 'p1';

% mn = 'ZYE_0089';
% td = '2024-07-09';
% en = 1;
% imec = 0;

serverRoot = expPath(mn, td, en);
%%
%plot PCA components and traces for blue and purple channels
corrPath = fullfile(serverRoot, 'corr', 'svdTemporalComponents_corr.npy');
if ~exist(corrPath, 'file')
    colors = {'blue', 'violet'};
    computeWidefieldTimestamps(serverRoot, colors); % preprocess video
    nSV = 200;
    [U, V, t, mimg] = hemoCorrect(serverRoot, nSV); % process hemodynamic correction
else
    nSV = 200;
    [U, V, t, mimg] = loadUVt(serverRoot, nSV);
end
if length(t) > size(V,2)
  t = t(1:size(V,2));
end
dV = [zeros(size(V,1),1) diff(V,[],2)];
ddV = [zeros(size(dV,1),1) diff(dV,[],2)];
tlFile = fullfile(serverRoot, 'blue\meanImage.npy');
meanImage = readNPY(tlFile);
%% 
figure;
imagesc(meanImage);
colormap(gray);
axis image;
point = [280,324];
% point = [214,353];
% point = [364,429];
% point = [376,448];
% point = [325,379];
hold on;
scatter(point(2),point(1),12,'r','filled');
%% realign the time of wf trace to ephys activity
syncTL = loadAlign(serverRoot, 'tl');
% syncProbe = loadAlign(serverRoot, [probeName '_imec' num2str(imec)]);
syncProbe = loadAlign(serverRoot, probeName);
if size(syncProbe,1)-size(syncTL,1) ~= 0
    % syncProbe = [syncProbe; zeros(size(syncTL,1)-size(syncProbe,1), 1)];
    syncTL = syncTL(1:size(syncProbe,1),1);
end
tAll1 = interp1(syncTL, syncProbe, t); 
%%
px1 = squeeze(U(point(1),point(2),:))'*dV;
px = px1/meanImage(point(1),point(2));
%%
figure;
plot(t,px(1:length(t)));
%%
clear dataMatrix2
epochSize = 20*35;
batchN = floor(size(px,2)/epochSize);
dataMatrix = zeros(batchN, epochSize);
for i = 1:batchN
    dataMatrix2(i,:) =  px(1,1+epochSize*(i-1):epochSize*i);
end
[freq1,psdx,psdx_mean] = fft_spectrum(dataMatrix2);
psdx_mean = mean(psdx,2);
psdx_std = std(psdx,[],2);
freq_value = [0.2,0.5,1,2,4,6,8];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean);
hs1g2 = figure('Renderer', 'painters', 'Position', [100 100 300 600]);
subplot(1,1,1);
plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
% shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),...
%     log10(psdx_std(2:freqN))./sqrt(size(psdx,2)),'lineProps','k');
hold on;
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (df/f^2)');
% ylim([-7,-3]);
%%
% epochT = [1100 1150];
% epochT = [1300 1350];
% epochT = [1000 1200];
epochT = [492 507];
epochT = interp1(syncTL, syncProbe, epochT);  
%%
indx1 = find(tAll1>epochT(1),1,'first');
indx2 = find(tAll1>epochT(2),1,'first');
figure;
plot(tAll1(indx1:indx2),px(indx1:indx2));

%%
ksRoot = fullfile(fileparts(getProbeFile(serverRoot, probeName)));
sp = loadKSdir(ksRoot);
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksRoot);
%%
figure;
incl1 = (spikeAmps>20);
spikeTimes2 = spikeTimes(incl1);
spikeDepths2 = spikeDepths(incl1);
plot(spikeTimes2, spikeDepths2,'.','color','k')
%%
nChansInFile = 385;
d = dir(fullfile(ksRoot,'*.bin'));
lfpFilename = fullfile(d.folder, d.name);
nSamps = d.bytes/2/nChansInFile;
mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
%%
ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
run(fullfile('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysWidefieldPhasemap', 'configure_file'))
chanMap1 = chanMapReorder(ops.chanMap);
%%
fs = 30000;
sampleEpoch = round(epochT)*fs;
chan = 90;
thisDat = double(mmf.Data.x(chan,sampleEpoch(1):sampleEpoch(2)));
fs_low = 20;
[b1, a1] = butter(3, fs_low/fs*2, 'low');
thisDat_filt = filtfilt(b1,a1,thisDat);
t1 = epochT(1):1/fs:epochT(2);
% downsample 300x to 100Hz
% t2 = epochT(1):1/fs_low:epochT(2);
% thisDat_filt2 = interp1(t1,thisDat_filt',t2);  
t2 = downsample(t1,300);
thisDat_filt2 = downsample(thisDat_filt,300);
%
figure;
incl1 = (spikeAmps>20 & spikeTimes>epochT(1) & spikeTimes<epochT(2));
spikeTimes2 = spikeTimes(incl1);
spikeDepths2 = spikeDepths(incl1);
plot(spikeTimes2, spikeDepths2,'.','color','k')
hold on;
plot(tAll1(indx1:indx2),px(indx1:indx2)*5000+1000,'color','r'); 
hold on;
time_max = max(spikeTimes2);
time_min = min(spikeTimes2);
binSize = 0.025;
[spikeHist,spikeT] = hist(spikeTimes2, time_min:binSize:time_max);
plot(spikeT,spikeHist+800,'color','g');
hold on;
plot(t2,thisDat_filt2/10+500,'m');
%%
epochSize = 2*fs_low;
batchN = floor(size(thisDat_filt2,2)/epochSize);
dataMatrix = zeros(batchN, epochSize);
for i = 1:batchN
    dataMatrix(i,:) =  thisDat_filt2(1,1+epochSize*(i-1):epochSize*i);
end
%%    
% trace is ntrace x timestamps
alpha_trace_mean = bsxfun(@minus, dataMatrix, mean(dataMatrix, 2)); % mean-subtract each SVD
N = size(alpha_trace_mean,2);
xdft = fft(alpha_trace_mean');
Fs = 100;
xdft1 = xdft(1:N/2+1,:);
nf = size(xdft1,1);
psdx = (1/(Fs*N)) * abs(xdft1).^2;
psdx(2:end-1,:) = 2*psdx(2:end-1,:);
%%
freq1 = 0:Fs/N:Fs/2;
psdx_mean = mean(psdx,2);
psdx_std = std(psdx,[],2);
freq_value = [0.2,0.5,1,2,4,6,8,10,16,30];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean);
hs1g2 = figure('Renderer', 'painters', 'Position', [100 100 300 600]);
subplot(1,1,1);
plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
% shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),...
%     log10(psdx_std(2:freqN))./sqrt(size(psdx,2)),'lineProps','k');
hold on;
% xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8','10','16','25'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (df/f^2)');
% ylim([-7,-3]);