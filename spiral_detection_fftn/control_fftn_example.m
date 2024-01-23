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
%% cortex surface outline
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%%
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
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
%%
params.downscale = 1;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
%% registration
fname1 = [mn '_' tdb '_' num2str(en) '_tform'];
regist_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration';
load(fullfile(regist_folder, fname1));
%%
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetection4\nullModelROI\roi';
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(roi_folder,[fname '_roi.mat']));
U1 = U./mimg;
U1 = U1(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2),1:50);
mimg1 = mimg(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2));
%%
tStart = 1681; tEnd = 1683; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(1:50,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U1,dV1,t,params,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%% 3d fft
[trace2d1_fft,traceAmp1_fft,tracePhase1_fft] = spiralPhaseMap_fftn(U1,dV1,t,params,rate);
trace2d1_fft = trace2d1_fft(:,:,1+35/rate:end-35/rate); 
tracePhase1_fft = tracePhase1_fft(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
a = 150; b = 512; c = 1; d = 512;
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
h1 = figure;
ax6 = subplot(1,1,1);
imagesc(squeeze(tracePhase1a(:,:,20)));
colormap(ax6,colorcet('C06'));
axis image; 
% axis off;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
%%
trace2 = zeros(512,512,size(trace2d1,3)); trace2(150:512,:,:) = trace2d1;
trace2_fft = zeros(512,512,size(trace2d1,3)); trace2_fft(150:512,:,:) = trace2d1_fft;
tracePhase2 = zeros(512,512,size(trace2d1,3)); tracePhase2(150:512,:,:) = tracePhase1;
tracePhase2_fft = zeros(512,512,size(trace2d1,3)); tracePhase2_fft(150:512,:,:) = tracePhase1_fft;
%%
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
dfolder = fullfile(folder,mn,td,num2str(en));
tformName = dir(fullfile(dfolder,'*tform.mat')).name;
load(fullfile(dfolder,tformName));
%%
trace3 = imwarp(trace2,tform,'OutputView',imref2d(size(projectedTemplate1)));
trace3_fft = imwarp(trace2_fft,tform,'OutputView',imref2d(size(projectedTemplate1)));
tracePhase3 = imwarp(tracePhase2,tform,'OutputView',imref2d(size(projectedTemplate1)));
tracePhase3_fft = imwarp(tracePhase2_fft,tform,'OutputView',imref2d(size(projectedTemplate1)));
%%
BW1 = logical(projectedAtlas1);
mask1 = zeros(512,512);
mask1(150:512,:) = 1;
mask2 = imwarp(mask1,tform,'OutputView',imref2d(size(projectedTemplate1)));
BW2 = logical(mask2);
BW = (BW1&BW2);
%%
dff = trace3*100; dff_fft = trace3_fft*100;
raw_min = min(dff_fft(:)); raw_max = max(dff_fft(:));
scale3 = 5;
scale2 = 1;
startframe = 5;
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 900]);
for i = 1:10
    ax1(i) = subplottight(4,11,i);
    im1 = imagesc(squeeze(dff(:,:,i+startframe)));
    set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    caxis([raw_min,raw_max]);
    colormap(ax1(i),parula);
    axis image; axis off;
    
    ax2(i) = subplottight(4,11,i+11);
    im2 = imagesc(squeeze(tracePhase3(:,:,i+startframe)));
    set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax2(i),colorcet('C06'));
    axis image; axis off;
    
    ax3(i) = subplottight(4,11,i+22);
    im3 = imagesc(squeeze(dff_fft(:,:,i+startframe)));
    set(im3, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    caxis([raw_min,raw_max]);
    colormap(ax3(i),parula);
    axis image; axis off;
    
    ax4(i) = subplottight(4,11,i+33);
    im4 = imagesc(squeeze(tracePhase3_fft(:,:,i+startframe)));
    set(im4, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath(1:3),st,atlas1,[],scale3);
    plotOutline(maskPath(4),st,atlas1,[],scale3);
    plotOutline(maskPath(5),st,atlas1,[],scale3);
    plotOutline(maskPath(6:11),st,atlas1,[],scale3);
    plotOutline(maskPath(1:11),st,atlas1,[],scale3);
    set(gca,'Ydir','reverse')
    colormap(ax4(i),colorcet('C06'));
    axis image; axis off;
end
i = 11;
ax1(11) = subplottight(4,11,i);
im1 = imagesc(squeeze(dff(:,:,i+startframe)));
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath(1:3),st,atlas1,[],scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
plotOutline(maskPath(1:11),st,atlas1,[],scale3);
set(gca,'Ydir','reverse')
caxis([raw_min,raw_max]);
colormap(ax1(i),parula);
axis image; axis off;
colorbar;

ax3(i) = subplottight(4,11,i+22);
im3 = imagesc(squeeze(dff_fft(:,:,i+startframe)));
set(im3, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath(1:3),st,atlas1,[],scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
plotOutline(maskPath(1:11),st,atlas1,[],scale3);
set(gca,'Ydir','reverse')
caxis([raw_min,raw_max]);
colormap(ax3(i),parula);
axis image; axis off;
colorbar
%%
print(h1, 'control_fftn_example', '-dpdf', '-bestfit', '-painters');