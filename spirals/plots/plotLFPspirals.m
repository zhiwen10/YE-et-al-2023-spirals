function [hs4bc,hs4de] = plotLFPspirals(data_folder, save_folder)
%% session info
mn = 'ZYE_0020';
td = '2021-03-26';
en = 1; imec = 0; probeName = 'p1';
% save lowpass filtered and dowsampled LFP data to ops.fproc1
load(fullfile(data_folder,'spirals','spirals_LFP','params.mat'));
ops.fproc1 = fullfile(data_folder,'spirals','spirals_LFP',...
    'filtered_lfp_imec0.dat');
ops.chanMap = fullfile(data_folder,'ephys','config_files',...
    'NPtype24_hStripe_botRow0_ref1.mat')
% total number of good channels that we will spike sort
%% check downsampled data with raw ephys data
epochT = [1215,1220];
ichan = 50;
lfpFs = 100;
%% plot ephys channel phasemap
Fs = 35;
[traceEphysMap1,phaseEphysMap1,ampEphysMap] = lfpPhasemap(ops, epochT, Fs);
chanMap1 = chanMapReorder(ops.chanMap);
%%
for j = 1:384
    [row(j),col(j)] = ind2sub(size(chanMap1),find(chanMap1 == j));
    phaseEphysMap(:,row(j),col(j)) = phaseEphysMap1(:,j);
    traceEphysMap(:,row(j),col(j)) = traceEphysMap1(:,j);
end
%%
traceEphysMap1 = traceEphysMap1*2.34/1000;                                 % multiply by amplifier gain
traceEphysMap = traceEphysMap*2.34/1000;                                   % multiply by amplifier gain
%% 
traceEphysSize = size(traceEphysMap1,1);
sizeN = size(phaseEphysMap,1);
immax = max(max(max(phaseEphysMap)));
immin = min(min(min(phaseEphysMap)));
maxTraceEphys = max(max(max(traceEphysMap)));
minTraceEphys = min(min(min(traceEphysMap)));
imH = []; imH2 = [];
sizeN = size(phaseEphysMap1,1);
immax = max(phaseEphysMap1(:)); immin = min(phaseEphysMap1(:));
maxRaw = max(traceEphysMap(:)); minRaw = min(traceEphysMap(:));
%%
% your data should be size [nFrames nChannels] where the channels are in
% the same order that the coordinates are in the channel map. so don't
% resize the data to be 8 x 48
chanMap = load(ops.chanMap);
caxPhase = [immin immax];
caxRaw = [minRaw/2 maxRaw/2];
colMapRaw = parula(100);
colMapPhase = colorcet('C06','N',100);
%%
hs4bc = figure;
hs4bc.Position = [100 100 800 500];
kk = 37;
ax1 = subplot(2,2,1);
set(ax1,'color','k');
psRaw = makeProbePlot(ax1, chanMap, 12);
colorSites(psRaw, traceEphysMap1(kk,:), colMapRaw, caxRaw);
axis off; 
ax3 = subplot(2,2,3);
set(ax3,'color','k');
psPhase = makeProbePlot(ax3, chanMap, 12);
colorSites(psPhase, phaseEphysMap1(kk,:), colMapPhase, caxPhase);
axis off;    

ax2 = subplot(2,2,2);
imH1 = imagesc(squeeze(traceEphysMap(kk,:,:)));
axis image; 
axis off;
colormap(ax2,parula); 
caxis(ax2,caxRaw);
set(ax2, 'ydir','normal');

ax4 = subplot(2,2,4);
imH1 = imagesc(squeeze(phaseEphysMap(kk,:,:)));
axis image; 
axis off;
colormap(ax4,colorcet('C06'));
set(ax4, 'ydir','normal');
caxis(ax4,caxPhase);
%%
print(hs4bc,fullfile(save_folder,'FigS4bc_LFP_4shank_example.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%%
hs4de = figure;
hs4de.Position = [100 100 800 500];
hs4de.Color = 'k';
count1 = 1;
for kk = 27:46
    axh2a(count1) = subplot(2,20,count1);
    imH1 = imagesc(squeeze(traceEphysMap(kk,:,:)));
    axis image; 
    axis off;
    colormap(axh2a(count1),parula); 
    caxis(caxRaw);
    set(axh2a(count1), 'ydir','normal')
    if count1 == 20
        colorbar;
    end

    axh2b(count1) = subplot(2,20,count1+20);
    imH1 = imagesc(squeeze(phaseEphysMap(kk,:,:)));
    axis image; 
    axis off;
    colormap(axh2b(count1),colorcet('C06'));
    set(axh2b(count1), 'ydir','normal')
    caxis(caxPhase);
    if count1 == 20
        colorbar;
    end
    
    count1 = count1+1;
end
%%
print(hs4de, fullfile(save_folder,'FigS4de_LFP_4shank_spirals.pdf'),...
    '-dpdf', '-bestfit', '-painters');