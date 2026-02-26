function hs3ab = plotExampleDataVsFft(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% session info
kk = 7;                                                                    % use ZYE_0012 as an example session
mn = T.MouseID{kk};
en = T.folder(kk);    
tda = T.date(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%%
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));       % load atlas transformation matrix tform
load(fullfile(data_folder,'spirals','fft_roi',[fname '_roi.mat']));          % read rectangle ROIouseID{kk};%% read svd components (U,V,t) from processed data folder
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];    
%% apply rectangle roi to the image
U1 = U./mimg;
U1 = U1(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2),1:50);
mimg1 = mimg(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2));
%%
params.downscale = 1;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
%%
tStart = 1681; tEnd = 1683; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(1:50,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%% 3d fft
[trace2d1_fft,traceAmp1_fft,tracePhase1_fft] = spiralPhaseMap_fftn_freq(U1,dV1,t,params,freq,rate);
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
%%
trace2 = zeros(512,512,size(trace2d1,3)); trace2(150:512,:,:) = trace2d1;
trace2_fft = zeros(512,512,size(trace2d1,3)); trace2_fft(150:512,:,:) = trace2d1_fft;
tracePhase2 = zeros(512,512,size(trace2d1,3)); tracePhase2(150:512,:,:) = tracePhase1;
tracePhase2_fft = zeros(512,512,size(trace2d1,3)); tracePhase2_fft(150:512,:,:) = tracePhase1_fft;
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
hs3ab = figure('Renderer', 'painters', 'Position', [100 100 900 900]);
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
colorbar;
%%
print(hs3ab, fullfile(save_folder,'FigS3ab_control_fftn_example.pdf'),...
    '-dpdf', '-bestfit', '-painters');