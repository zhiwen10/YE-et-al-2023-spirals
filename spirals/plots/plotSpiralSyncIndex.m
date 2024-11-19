function [hs8e, hs8ad,hs8f] = plotSpiralSyncIndex(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
kk = 5;                                                                    % ZYE_0060
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root); 
dV = [zeros(size(V,1),1) diff(V,[],2)];
%%
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform;
load(fullfile(data_folder,'spirals','spirals_index',[fname '_motion_energy.mat']));
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
freq = [2,8];
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
BW = logical(projectedAtlas1);
BW = BW(1:params.downscale:end,1:params.downscale:end);
%%
load(fullfile(data_folder,'spirals','spirals_grouping',...
    [fname '_spirals_group_fftn.mat']));                                   % load spiral centers (>40 pixels radius)
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform
spiral_duration = cellfun(@(x) size(x,1), archiveCell);
indx2 = (spiral_duration>=15);
groupedCells = archiveCell(indx2);
%%
frameStart = 72700; frameEnd = frameStart+351;
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter
traceAmp1 = traceAmp1(:,:,1+35/rate:end-35/rate)./mimgtransformed; % reduce 2*35 frames after filter data 
%%
trace1 = squeeze(Utransformed(100,100,1:50))'*dV(1:50,:);
trace1 = trace1./mimgtransformed(100,100);
%%
hemi_select = 72:143;
cols = size(tracePhase1,1);
rows= size(tracePhase1,2);
xs = ceil(linspace(-rows/4,rows/4,ceil(rows/2)));
ys = ceil(linspace(-cols/2,cols/2,cols));
[Xs,Ys] = meshgrid(xs,ys);
anglein = atan2(Ys,Xs);  % expected angle for centered sprial
BW_half = BW(:,hemi_select);
%%
clear amp_mean1 index_unity1 sync1 spirality1
for frame = 1:size(tracePhase1,3)
    phase_hemi = squeeze(tracePhase1(:,hemi_select,frame));
    phase_hemi(not(BW_half)) = nan;
    pixel_count = sum(BW_half(:));
    sync1(frame,1)  = abs(sum(exp(j*phase_hemi(:)),"omitnan"))./pixel_count;
    phase_diff = phase_hemi-anglein;
    spirality1(frame,1) = abs(sum(exp(j*phase_diff(:)),"omitnan"))./pixel_count;
    amp_hemi = squeeze(traceAmp1(:,hemi_select,frame));
    amp_hemi(not(BW_half)) = nan;
    amp_hemi = amp_hemi(:);
    amp_mean1(frame,1) = mean(amp_hemi,"omitnan");
end
index_unity1 = sqrt(sync1.^2+spirality1.^2);
motion_energy = image_energy2(frameStart:frameEnd);
%% plot example time series 
frame1 = 42; frame2 = 70;
t2 = size(amp_mean1,1);
hs8e = figure;
subplot(6,1,1);
plot(1/35:1/35:t2/35,trace1(frameStart:frameEnd));
xline(frame1/35,'--');
xline(frame2/35,'--');
subplot(6,1,2);
plot(1/35:1/35:t2/35,amp_mean1);
ylabel('3-6Hz amp');
xline(frame1/35,'--');
xline(frame2/35,'--');

subplot(6,1,3);
plot(1/35:1/35:t2/35,motion_energy);
ylabel('motion energy');
xline(frame1/35,'--');
xline(frame2/35,'--');

subplot(6,1,4);
plot(1/35:1/35:t2/35,sync1);
ylim([0,1]);
ylabel('sync index');
xline(frame1/35,'--');
xline(frame2/35,'--');

subplot(6,1,5);
plot(1/35:1/35:t2/35,spirality1);
ylim([0,1]);
ylabel('spirarity index');
xline(frame1/35,'--');
xline(frame2/35,'--');

subplot(6,1,6);
plot(1/35:1/35:t2/35,index_unity1);
ylim([0,1]);
ylabel('sum index');
xline(frame1/35,'--');
xline(frame2/35,'--');

print(hs8e, fullfile(save_folder,'Figs8e_example_time_series'),...
    '-dpdf', '-bestfit', '-painters');
%% plot example frames for index
frame1 = 42; frame2 = 70;
hs8ad = figure('Renderer', 'painters', 'Position', [100 100 900 360]);
ax1 = subplot(2,6,3);
phasemap1 = squeeze(tracePhase1(:,hemi_select,frame1));
phasemap1(not(BW_half)) = nan;
im_phase = imagesc(phasemap1);
colormap(ax1,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax2 = subplot(2,6,4);
im_phase = imagesc(anglein);
colormap(ax2,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax3 = subplot(2,6,5);
angle_diff1 = phasemap1-anglein;
angle_diff1 = wrapToPi(angle_diff1);
im_phase = imagesc(angle_diff1);
colormap(ax3,colorcet('C06'));
% caxis([-pi,pi]);
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;

ax4 = subplot(2,6,6);
angle_diff2 = angle_diff1(~isnan(angle_diff1));
angle_diff2_cV = exp(j*angle_diff2);
cmean = circ_mean(angle_diff2);
cvariance = circ_r(angle_diff2);
cmean_cV = exp(j*cmean)*cvariance;
x1 = real(angle_diff2_cV);
y1 = imag(angle_diff2_cV);
p = randperm(numel(x1),500);
x2 = x1(p); y2 = y1(p);
x = zeros(size(x2));
y = zeros(size(y2));
quiver(x,y,x2,y2,'k');
hold on;
quiver(0,0,real(cmean_cV),imag(cmean_cV),'r');
axis image;
xlim([-1,1]);
ylim([-1,1]);
text(-1,0.9,['spiralIndex = ', num2str(round(cvariance*10)/10)]);

ax5 = subplot(2,6,1);
im_phasemap = squeeze(tracePhase1(:,hemi_select,frame1));
im_phasemap(not(BW_half)) = nan;
im_phase = imagesc(im_phasemap);
colormap(ax5,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax6 = subplot(2,6,2);
angle_diff3 = im_phasemap(~isnan(im_phasemap));
angle_diff3_cV = exp(j*angle_diff3);
cmean1 = circ_mean(angle_diff3);
cvariance1 = circ_r(angle_diff3);
cmean_cV1 = exp(j*cmean1)*cvariance1;
x1 = real(angle_diff3_cV);
y1 = imag(angle_diff3_cV);
p = randperm(numel(x1),500);
x2 = x1(p); y2 = y1(p);
x = zeros(size(x2));
y = zeros(size(y2));
quiver(x,y,x2,y2,'k');
hold on;
quiver(0,0,real(cmean_cV1),imag(cmean_cV1),'r');
axis image;
xlim([-1,1]);
ylim([-1,1]);
text(-1,0.9,['syncIndex = ', num2str(round(cvariance1*10)/10)]);
%
ax1 = subplot(2,6,9);
phasemap1 = squeeze(tracePhase1(:,hemi_select,frame2));
phasemap1(not(BW_half)) = nan;
im_phase = imagesc(phasemap1);
colormap(ax1,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax2 = subplot(2,6,10);
im_phase = imagesc(anglein);
colormap(ax2,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax3 = subplot(2,6,11);
angle_diff1 = phasemap1-anglein;
angle_diff1 = wrapToPi(angle_diff1);
im_phase = imagesc(angle_diff1);
colormap(ax3,colorcet('C06'));
% caxis([-pi,pi]);
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;

ax4 = subplot(2,6,12);
angle_diff2 = angle_diff1(~isnan(angle_diff1));
angle_diff2_cV = exp(j*angle_diff2);
cmean = circ_mean(angle_diff2);
cvariance = circ_r(angle_diff2);
cmean_cV = exp(j*cmean)*cvariance;
x1 = real(angle_diff2_cV);
y1 = imag(angle_diff2_cV);
p = randperm(numel(x1),500);
x2 = x1(p); y2 = y1(p);
x = zeros(size(x2));
y = zeros(size(y2));
quiver(x,y,x2,y2,'k');
hold on;
quiver(0,0,real(cmean_cV),imag(cmean_cV),'r');
axis image;
xlim([-1,1]);
ylim([-1,1]);
text(-1,0.9,['spiralIndex = ', num2str(round(cvariance*10)/10)]);

ax5 = subplot(2,6,7);
im_phasemap = squeeze(tracePhase1(:,hemi_select,frame2));
im_phasemap(not(BW_half)) = nan;
im_phase = imagesc(im_phasemap);
colormap(ax5,colorcet('C06'));
set(im_phase, 'AlphaData', BW_half, 'AlphaDataMapping', 'scaled');
axis off; axis image;
caxis([-pi,pi]);

ax6 = subplot(2,6,8);
angle_diff3 = im_phasemap(~isnan(im_phasemap));
angle_diff3_cV = exp(j*angle_diff3);
cmean1 = circ_mean(angle_diff3);
cvariance1 = circ_r(angle_diff3);
cmean_cV1 = exp(j*cmean1)*cvariance1;
x1 = real(angle_diff3_cV);
y1 = imag(angle_diff3_cV);
p = randperm(numel(x1),500);
x2 = x1(p); y2 = y1(p);
x = zeros(size(x2));
y = zeros(size(y2));
quiver(x,y,x2,y2,'k');
hold on;
quiver(0,0,real(cmean_cV1),imag(cmean_cV1),'r');
axis image;
xlim([-1,1]);
ylim([-1,1]);
text(-1,0.9,['syncIndex = ', num2str(round(cvariance1*10)/10)]);
%
print(hs8ad, fullfile(save_folder,'Figs8ad_example_index.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%% plot sync and spirality relationship
frameStart = 72700; frameEnd = frameStart+3501;
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter
traceAmp1 = traceAmp1(:,:,1+35/rate:end-35/rate)./mimgtransformed; % reduce 2*35 frames after filter data 
clear amp_mean1 index_unity1 sync1 spirality1
for frame = 1:size(tracePhase1,3)
    phase_hemi = squeeze(tracePhase1(:,hemi_select,frame));
    phase_hemi(not(BW_half)) = nan;
    pixel_count = sum(BW_half(:));
    sync1(frame,1)  = abs(sum(exp(j*phase_hemi(:)),"omitnan"))./pixel_count;
    phase_diff = phase_hemi-anglein;
    spirality1(frame,1) = abs(sum(exp(j*phase_diff(:)),"omitnan"))./pixel_count;
    amp_hemi = squeeze(traceAmp1(:,hemi_select,frame));
    amp_hemi(not(BW_half)) = nan;
    amp_hemi = amp_hemi(:);
    amp_mean1(frame,1) = mean(amp_hemi,"omitnan");
end
index_unity1 = sqrt(sync1.^2+spirality1.^2);
motion_energy = image_energy2(frameStart:frameEnd);
%
color2 = flipud(cbrewer2('div','RdBu',32));
% color2 = cbrewer2('seq','Greys',32);
hs8f = figure;
scatter(sync1,spirality1,16,log10(amp_mean1),'filled');
colormap(color2)
hold on;
r = 1;
th = 0:pi/100:pi/2;
xunit = r * cos(th);
yunit = r * sin(th);
plot(xunit, yunit,'--k');
xlabel('synchrony');
ylabel('spirality');
axis image
xlim([0,1]);
ylim([0,1]);
colorbar
print(hs8f, fullfile(save_folder,'Figs8f_sync_spirality_amp'),...
    '-dpdf', '-bestfit', '-painters');