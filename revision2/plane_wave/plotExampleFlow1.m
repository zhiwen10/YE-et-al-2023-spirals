function plotExampleFlow1(data_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% load data from an example an session
mn = 'ZYE_0012';
td = '2020-10-16';
tdb = datestr(td,'yyyymmdd');
en = 5;
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'tables',[fname '_tform.mat']));                 % load atlas transformation matrix tform;
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
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
freq = [2,8];
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
scale3 = 5/8;
color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
frame = 30;
phase_min = -pi; phase_max = pi;
%%
h1 = figure;
ax1 = subplot(1,1,1);
imagesc(tracePhase1(:,:,frame));
colormap(ax1,colorcet('C06'));
axis image;
%%
lineColor = 'k'; lineColor1 = 'w';
hemi = [];

h1 = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
framea = tracePhase1(:,:,frame);
frameb = tracePhase1(:,:,frame+1);

ax1 = subplottight(1,4,1);
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
% v = [70 50; 130 50;130,120;70 120];
% v = [100 60; 130 60;130,90;100 90];
v = [70 60; 130 60;130,120;70 120];
f = [1,2,3,4];
patch('Faces',f,'Vertices',v,'EdgeColor','w','FaceColor','none','LineWidth',2);
% v2 = [70 90; 100 90;100,120;70 120];
% patch('Faces',f,'Vertices',v2,'EdgeColor','k','FaceColor','none','LineWidth',2);

ax2 = subplottight(1,4,2);
im_phase = imagesc(frameb);
colormap(ax2,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
patch('Faces',f,'Vertices',v,'EdgeColor','w','FaceColor','none','LineWidth',2);
% patch('Faces',f,'Vertices',v2,'EdgeColor','k','FaceColor','none','LineWidth',2);

ax3 = subplottight(1,4,3);
frame_ab(1,:,:) = framea;
frame_ab(2,:,:) = frameb;
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(frame_ab,useGPU);
vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
skip = 4;
zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW2) = nan; vyRaw2(~BW2) = nan;

im_phase = imagesc(framea);
colormap(ax3,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
patch('Faces',f,'Vertices',v,'EdgeColor','w','FaceColor','none','LineWidth',2);

ax6 = subplottight(1,4,4);
im_phase = imagesc(framea);
colormap(ax6,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
% xlim([70,100]);
% ylim([90,120]);
xlim([70,130]);
ylim([60,120]);
ax6.Position(2) = ax6.Position(2)+0.05;
ax6.Position(3) = ax6.Position(3)-0.05;
ax6.Position(4) = ax6.Position(4)-0.05;
%%
print(h1,'exampleFlow1.pdf','-dpdf', '-bestfit', '-painters');