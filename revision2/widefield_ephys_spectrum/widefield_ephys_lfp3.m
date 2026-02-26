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
ops.chanMap = 'C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master\configFiles\NPtype24_hStripe_botRow0_ref1.mat';
run(fullfile('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\ephysWidefieldPhasemap', 'configure_file'))
chanMap1 = chanMapReorder(ops.chanMap);
%%
epochT = [495 505];
point(4,:) = [217,312]; point(3,:) = [218,326]; point(2,:) = [215,340]; point(1,:) = [214,353];
% ops.fproc1 = 'E:\ephys\ZYE_0016\2021-02-25\lfp_filtered_downSample_all.dat';
tt1 = epochT(1); tt2 = epochT(2);
%% realign the time of wf trace to ephys activity
syncTL = loadAlign(serverRoot, 'tl');
% syncProbe = loadAlign(serverRoot, [probeName '_imec' num2str(imec)]);
syncProbe = loadAlign(serverRoot, probeName);
if size(syncProbe,1)-size(syncTL,1) ~= 0
    % syncProbe = [syncProbe; zeros(size(syncTL,1)-size(syncProbe,1), 1)];
    syncTL = syncTL(1:size(syncProbe,1),1);
end
tAll1 = interp1(syncTL, syncProbe, t);  
%% filter wf trace at 2-8Hz with filtfilt
for i  = 1:4
    px_mean = mimg(point(i,1),point(i,2));
    px1 = squeeze(U(point(i,1),point(i,2),:))'*V;
    px(i,:) = px1/px_mean;
end
Fs = 35;
px = double(px);
% [f1,f2] = butter(2, [2 8]/(Fs/2), 'bandpass');
% pxTemp = filtfilt(f1,f2,px');
% px = pxTemp';
%%
% load spikes
ksRoot = fullfile(fileparts(getProbeFile(serverRoot, probeName)));
sp = loadKSdir(ksRoot);
[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(ksRoot);
%%
% shank sites
if contains(ops.chanMap,'NPtype24_hStripe')
    shanks{1} = find(sp.xcoords<-200);
    shanks{2} = find(sp.xcoords<100 & sp.xcoords>-200);
    shanks{3} = find(sp.xcoords<300 & sp.xcoords>100);
    shanks{4} = find(sp.xcoords>300);
elseif contains(ops.chanMap,'NPtype24_doubleLengthStripe')
    shanks{1} = find(sp.xcoords<100);
    shanks{2} = find(sp.xcoords<300 & sp.xcoords>200);
    shanks{3} = find(sp.xcoords<600 & sp.xcoords>400);
    shanks{4} = find(sp.xcoords>700);
end
%% MUA histogram and plot for 4 shanks
tAll2 = tAll1(not(isnan(tAll1)));
px2 = px(:,not(isnan(tAll1)));
twf = tt1:1/35:tt2;
pxTemp = interp1(tAll2,px2',twf);
pxTemp = pxTemp';
binSize = 0.025;
color1 = {'k','r','g','c'};
for shank = 1:4
    clear spikeT spikeHist 
    incl1 = (spikeAmps>20 & ismember(spikeSites,...
            shanks{shank}));%& spikeTimes>(3000/35) & spikeTimes<(3500/35));
    s1.spikeTimes = spikeTimes(incl1);
    s1.spikeDepths = spikeDepths(incl1);
    s2{shank} = s1;
    [spikeHist,spikeT] = spike_histogram(s2{shank}.spikeTimes,binSize);
    sIndex  = find(spikeT>=tt1 & spikeT <=tt2);
    spikeT = spikeT(sIndex);
    spikeHist = spikeHist(sIndex);
end
%% load raw ephys data
NchanTOT = ops.NchanTOT;
fs = ops.fs;
tTotal = tt2-tt1;
NTbuff = tTotal*fs;
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
StartSample = floor(tt1*fs);
offset1 = 2*NchanTOT*StartSample; % number of samples to start reading at.
fseek(fid, offset1, 'bof'); % fseek to batch start in raw file
buffa = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
fclose(fid)
Map = ops.chanMap;
% [chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(Map); % function to load channel map file
buffa  = double(buffa(ops.chanMap1,:)); % subsample only good channels
buffaMean = mean(buffa,1);
buffa = buffa-buffaMean;
chanMap1 = chanMapReorder(Map);
probeTips = chanMap1(3,[1,2,3,4]);
% probeTips = chanMap1(20,:);
ta = linspace(tt1,tt2,size(buffa,2));
%
fs_low = 100;
[b1, a1] = butter(3, fs_low/fs*2, 'low');
% thisDat_filt = filtfilt(b1,a1,thisDat);
thisDat_filt = filtfilt(b1,a1,buffa);
thisDat_filt2 = thisDat_filt(probeTips(shank),:);
t1 = epochT(1):1/fs:epochT(2);
t2 = epochT(1):1/fs_low:epochT(2);
thisDat_filt2 = interp1(t1(1:end-1),thisDat_filt2,t2);  
tb = tt1:1/100:tt2;
f1 = figure;
f1.Position = [100 100 1500 900];
%
shank = 1;
clear spikeT spikeHist incl1
% plot MUA activity
incl1 = (spikeAmps>20 & ismember(spikeSites,...
        shanks{shank})); %spikeTimes>(3000/35) & spikeTimes<(3500/35));
s1.spikeTimes = spikeTimes(incl1);
s1.spikeDepths = spikeDepths(incl1);
[spikeHist,spikeT] = spike_histogram(s1.spikeTimes,binSize);
sIndex  = find(spikeT>=tt1 & spikeT <=tt2);
spikeT = spikeT(sIndex);
spikeHist = spikeHist(sIndex);
%
% plot widefield trace at probe site
pxTrace = pxTemp(shank,:);
% plot LFP data
% buffDown = buffTemp4(probeTips(shank),tt1*100+1:floor(tt2*100));
% buffDown = thisDat_filt2(tt1*100+1:floor(tt2*100));
subplot(4,2,[1,3,5,7])
imagesc(mimg);
colormap(gray)
hold on;
scatter(point(shank,2),point(shank,1),'r');
axis image
axis off
subplot(4,2,2)
plot(twf, pxTrace*100,'c','lineWidth',1); 
legend('cortexWF')
xticks(tt1:2:tt2)
xlabel('Time (s)'); ylabel('df/f %');
xlim([tt1 tt2]);

subplot(4,2,4)
% spikeHist = ones(numel(spikeT),1);
plot(spikeT,spikeHist ,'m','lineWidth',1)
xlim([tt1 tt2]);

subplot(4,2,6)
% plot(spikeT,spikeHist/100+2 ,'m','lineWidth',1)
%plot raw data at the 4 shank sites
hold on; plot(ta,buffa(probeTips(shank),:)/500,'k');
hold on; 
plot(tb,thisDat_filt2/500,'r','lineWidth',1);
% legend('MUA','Raw','LFP')
xticks(tt1:2:tt2)
xlabel('Time (s)'); ylabel('uV');
xlim([tt1 tt2]);
ylim([-2,2]);

subplot(4,2,8)
plot(twf, pxTrace*20000+800,'c','lineWidth',1); 
incl0 = (spikeAmps>20 & ismember(spikeSites,...
        shanks{shank}) & spikeTimes>(tt1) & spikeTimes<(tt2));
spikeTimes1 = spikeTimes(incl0);
spikeDepths1 = spikeDepths(incl0);
hold on;
scatter(spikeTimes1, spikeDepths1,1,'k','filled')
xticks(tt1:2:tt2)
yticks(0:200:1000)
xlabel('Time (s)');
ylabel('ProbeDepth (um)');
xlim([tt1 tt2]);
%%
print(f1, 'widefield&ephys_ZYE16', '-dpdf', '-bestfit', '-painters');
%% get lfp long
% epochT = [450 500];
NchanTOT = ops.NchanTOT;
fs = ops.fs;
tTotal = epochT(2)-epochT(1);
NTbuff = tTotal*fs;
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
StartSample = floor(tt1*fs);
offset1 = 2*NchanTOT*StartSample; % number of samples to start reading at.
fseek(fid, offset1, 'bof'); % fseek to batch start in raw file
buffa = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
fclose(fid)
Map = ops.chanMap;
% [chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(Map); % function to load channel map file
buffa  = double(buffa(ops.chanMap1,:)); % subsample only good channels
buffaMean = mean(buffa,1);
buffa = buffa-buffaMean;
chanMap1 = chanMapReorder(Map);
probeTips = chanMap1(3,[1,2,3,4]);
% probeTips = chanMap1(20,:);
ta = linspace(epochT(1),epochT(2),size(buffa,2));
fs_low = 100;
[b1, a1] = butter(3, fs_low/fs*2, 'low');
% thisDat_filt = filtfilt(b1,a1,thisDat);
thisDat_filt = filtfilt(b1,a1,buffa);
thisDat_filt2 = thisDat_filt(probeTips(shank),:);
t1 = epochT(1):1/fs:epochT(2);
t2 = epochT(1):1/fs_low:epochT(2);
thisDat_filt2 = interp1(t1(1:end-1),thisDat_filt2,t2);  
%% get lfp long-2
% epochT = [400 600];
% nChansInFile = 385;
% d = dir(fullfile(ksRoot,'*.bin'));
% lfpFilename = fullfile(d.folder, d.name);
% nSamps = d.bytes/2/nChansInFile;
% mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});
% fs = 30000;
% sampleEpoch = round(epochT)*fs;
% chan = probeTips(1);
% thisDat = double(mmf.Data.x(chan,sampleEpoch(1):sampleEpoch(2)));
% fs_low = 100;
% [b1, a1] = butter(3, fs_low/fs*2, 'low');
% thisDat_filt = filtfilt(b1,a1,thisDat);
% t1 = epochT(1):1/fs:epochT(2);
% % downsample 300x to 100Hz
% t2 = downsample(t1,300);
% thisDat_filt2 = downsample(thisDat_filt,300);
%%
f2 = figure('Renderer', 'painters', 'Position', [100 100 600 400]);
subplot(1,2,1);
epochSize = 2*fs_low;
batchN = floor(size(thisDat_filt2,2)/epochSize);
dataMatrix = zeros(batchN, epochSize);
for i = 1:batchN
    dataMatrix(i,:) =  thisDat_filt2(1,1+epochSize*(i-1):epochSize*i);
end   
alpha_trace_mean = bsxfun(@minus, dataMatrix, mean(dataMatrix, 2)); % mean-subtract each SVD
N = size(alpha_trace_mean,2);
xdft = fft(alpha_trace_mean');
Fs = 100;
xdft1 = xdft(1:N/2+1,:);
nf = size(xdft1,1);
psdx = (1/(Fs*N)) * abs(xdft1).^2;
psdx(2:end-1,:) = 2*psdx(2:end-1,:);
freq1 = 0:Fs/N:Fs/2;
psdx_mean = mean(psdx,2);
psdx_std = std(psdx,[],2);
freq_value = [0.2,0.5,1,2,4,6,8,10,16,32];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean);
plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8','10','16','32'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (\muV)^2)');
% ylim([-7,-3]);

twf2 = epochT(1):1/35:epochT(2);
pxTemp = interp1(tAll2,px2',twf2);
pxTrace = pxTemp(:,1);
clear dataMatrix2
% pxTrace = px2(1,:);
epochSize = 2*35;
batchN = floor(length(pxTrace)/epochSize);
dataMatrix = zeros(batchN, epochSize);
for i = 1:batchN
    dataMatrix2(i,:) =  pxTrace(1+epochSize*(i-1):epochSize*i);
end
[freq1,psdx,psdx_mean2] = fft_spectrum(dataMatrix2);
psdx_mean2 = mean(psdx,2);
psdx_std = std(psdx,[],2);
freq_value = [0.2,0.5,1,2,4,6,8,10,16,32];
log_freq_value = log10(freq_value);
freqN = length(psdx_mean2);
subplot(1,2,2);
plot(log10(freq1(2:freqN)),log10(psdx_mean2(2:freqN)),'color','k');
% shadedErrorBar(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),...
%     log10(psdx_std(2:freqN))./sqrt(size(psdx,2)),'lineProps','k');
hold on;
xlim([log_freq_value(1),log_freq_value(end)]);
xticks(log_freq_value(1:end));
xticklabels({'0.2','0.5','1','2','4','6','8','10','16','32'});
xlabel('log10(Frequency)');
ylabel('log10(Power) (dF/F)^2)');
print(f2, 'widefield&ephys_ZYE16_power_spectrum_2', '-dpdf', '-bestfit', '-painters');
%% spike_histogram
function [n,x] = spike_histogram(spikeT,binSize)
time_max = max(spikeT);
time_min = min(spikeT);
[n,x] = hist(spikeT, time_min:binSize:time_max);
end