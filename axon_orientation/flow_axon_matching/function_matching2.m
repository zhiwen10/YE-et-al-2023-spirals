%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\nrrdread'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
mn = 'ZYE_0012';
td = '2020-10-16';
tdb = datestr(td,'yyyymmdd');
en = 5;
serverRoot = expPath(mn, td, en);
wf_svd
%% get atlas mask and outlines
% load coords for atlas
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
% load projectedAtlas and projectedTemplate
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%% registration
fname = [mn '_' tdb '_' num2str(en) '_tform'];
regist_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration';
load(fullfile(regist_folder, fname));
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
BW = logical(projectedAtlas1);
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
%%
[row,col] = find(BW);
brain_index = [col,row];
%%
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(Utransformed,dV1,t,params,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
tracePhase_current = tracePhase1(:,:,20:32);
tracePhase_current = permute(tracePhase_current,[3 1 2]);
[vxRaw,vyRaw] = HS_flowfield(tracePhase_current,useGPU);
%%
figure;
for i = 1:12
    subplot(3,4,i);
    im = imagesc(squeeze(tracePhase_current(i,:,:)));
    colormap(colorcet('C06'));
    axis image; axis off;
    set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
    hold on;
    vxRaw1 = squeeze(vxRaw(i,:,:)); vyRaw1 = squeeze(vyRaw(i,:,:));
    vxRaw2 = nan(size(vxRaw1));
    vyRaw2 = nan(size(vyRaw1));
    skip = 3; zoom_scale = 2;
    vxRaw2(1:skip:end,1:skip:end) = vxRaw1(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw1(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
    hold on;
    imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
end
   
%%
vxRaw_mean = squeeze(mean(vxRaw,1)); vyRaw_mean= squeeze(mean(vyRaw,1));
figure;
ax1 = subplot(1,1,1);
im = imagesc(squeeze(tracePhase_current(1,:,:)));
colormap(ax1,colorcet('C06'));
axis image; axis off;
hold on;
vxRaw2 = nan(size(vxRaw_mean));
vyRaw2 = nan(size(vyRaw_mean));
skip = 3; zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
