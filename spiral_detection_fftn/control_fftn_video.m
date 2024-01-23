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
%%
params.downscale = 1;
params.lowpass = 0;
params.gsmooth = 0;
rate = 0.1;
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
tStart = 1681; tEnd = 1685; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(1:50,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U1,dV1,t,params,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%% 3d fft
[trace2d1_fft,traceAmp1_fft,tracePhase1_fft] = spiralPhaseMap_fftn3(U1,dV1,t,params,rate);
trace2d1_fft = trace2d1_fft(:,:,1+35/rate:end-35/rate); 
tracePhase1_fft = tracePhase1_fft(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
dff = trace2d1*100;
dff_fft = trace2d1_fft*100;
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
v = VideoWriter(fname);
v.FrameRate = 30;
open(v);
raw_min = min(dff_fft(:)); raw_max = max(dff_fft(:));
phase_min = -3.14; phase_max = 3.14;
ha = figure('Renderer', 'painters', 'Position', [100 100 1000 800]);
set(gcf,'color','black')
ha.MenuBar = 'none';
ha.ToolBar = 'none';
ha.Resize = 'off';
ha.Clipping = 'off';

ax1 = subplot(2,2,1);
im_raw = imagesc(dff(:,:,1));
colormap(ax1, parula)
axis image; axis off;
% set(im_raw, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
cb1 = colorbar; 
cb1.FontSize = 16;
cb1.Color = 'w';
pos1 = get(cb1,'Position');
set(cb1,'position',[pos1(1)+0.08 pos1(2)+0.1 pos1(3) pos1(4)/2]);
caxis([raw_min,raw_max]);
text_ax1 = text(0,-20,[num2str(round(qt(1)*10)/10) ' s'],'fontSize',16,'color','w');
text1_ax1 = text(150,-50,'Raw trace (df/f, %)','fontSize',16,'color','w');

ax2 = subplot(2,2,2);
im_phase = imagesc(tracePhase1(:,:,1));
% colormap(ax4,hsv);
colormap(ax2,colorcet('C06'));
axis image; axis off;
% set(im_phase, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
cb2 = colorbar('Ticks',[-3.14, 0, 3.14],...
         'TickLabels',{'-pi','0','pi'});  
cb2.FontSize = 16;
cb2.Color = 'w';
pos2 = get(cb2,'Position');
set(cb2,'position',[pos2(1)+0.08 pos2(2)+0.1 pos2(3) pos2(4)/2]);
caxis([phase_min,phase_max]);
text2_ax1 = text(150,-50,'Raw phase','fontSize',16,'color','w');
%
ax3 = subplot(2,2,3);
im_raw_fft = imagesc(dff_fft(:,:,1));
colormap(ax3, parula)
axis image; axis off;
% set(im_raw, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
cb1 = colorbar; 
cb1.FontSize = 16;
cb1.Color = 'w';
pos1 = get(cb1,'Position');
set(cb1,'position',[pos1(1)+0.08 pos1(2)+0.1 pos1(3) pos1(4)/2]);
caxis([raw_min,raw_max]);
text3_ax1 = text(150,-50,'3D-FFT trace (df/f, %)','fontSize',16,'color','w');

ax4 = subplot(2,2,4);
im_phase_fft = imagesc(tracePhase1_fft(:,:,1));
% colormap(ax4,hsv);
colormap(ax4,colorcet('C06'));
axis image; axis off;
% set(im_phase, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
cb2 = colorbar('Ticks',[-3.14, 0, 3.14],...
         'TickLabels',{'-pi','0','pi'});  
cb2.FontSize = 16;
cb2.Color = 'w';
pos2 = get(cb2,'Position');
set(cb2,'position',[pos2(1)+0.08 pos2(2)+0.1 pos2(3) pos2(4)/2]);
caxis([phase_min,phase_max]);
text4_ax1 = text(150,-50,'3D-FFT phase','fontSize',16,'color','w');

start_frame_new = 1;
end_frame_new = numel(qt);
for i = start_frame_new:end_frame_new 
    B1t = squeeze(dff(:,:,i));
    set(im_raw, 'CData', B1t);   
    B1t3 = squeeze(tracePhase1(:,:,i));
    set(im_phase, 'CData', B1t3); 
    B1t_fft = squeeze(dff_fft(:,:,i));
    set(im_raw_fft, 'CData', B1t_fft);   
    B1t3_fft = squeeze(tracePhase1_fft(:,:,i));
    set(im_phase_fft, 'CData', B1t3_fft); 
    
    text_ax1.String = [num2str(round(qt(i)*10)/10) ' s'];
    thisFrame = getframe(ha); 
    writeVideo(v, thisFrame);
end
close(v);   
%%
rate1 = 1;
[trace2d1a,traceAmp1a,tracePhase1a] = spiralPhaseMap4(U(:,:,1:50),dV1,t,params,rate1);
trace2d1a = trace2d1a(:,:,1+35/rate1:end-35/rate1); 
tracePhase1a = tracePhase1a(:,:,1+35/rate1:end-35/rate1); % reduce 2*35 frames after filter data 
%%
a = 150; b = 512; c = 1; d = 512;
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
h1 = figure;
ax6 = subplot(1,1,1)
imagesc(squeeze(tracePhase1a(:,:,20)));
colormap(ax6,colorcet('C06'));
axis image; 
% axis off;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
%%
print(h1, 'ZYE12_FFT_ROI', '-dpdf', '-bestfit', '-painters');
%%
pixel(1,:) = [310,300]; % VISp
pixel(2,:) = [150,350]; % RSP
pixel(3,:) = [100,420]; % SSp-ul
%%
trace = zeros(3,size(trace2d1,3));
for i = 1:3
    % trace(i,:) = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'*dV1;
    trace(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
end
%%
trace_fft = zeros(3,size(trace2d1_fft,3));
for i = 1:3
    % trace(i,:) = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'*dV1;
    trace_fft(i,:) = trace2d1_fft(pixel(i,1),pixel(i,2),:);
end
%%
color1 = cbrewer2('qual','Set1',3,'spline');
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
figure;
ax1 = subplot(1,2,1);
for k = 1:3
    plot(ax1,qt,trace(k,:)*100+2*k,'color','k','lineWidth',1.5);
    hold on; 
end
xlim([t1(1),t1(end)]);
% set(ax1,'YTickLabel',[]);
% text(qt(1)-1,2,'VISp','color',color1(1,:),'fontSize',16,'Interpreter','none');
% text(qt(1)-1,4,'SSp-ul','color',color1(2,:),'fontSize',16,'Interpreter','none');
% text(qt(1)-1,6,'SSp-m','color',color1(3,:),'fontSize',16,'Interpreter','none');
xlabel('Time (s)','fontSize',16);
ax2 = subplot(1,2,2);
for k = 1:3
    plot(ax2,qt,trace_fft(k,:)*100+2*k,'color','k','lineWidth',1.5);
    hold on; 
end
xlim([t1(1),t1(end)]);
% set(ax1,'YTickLabel',[]);
% text(qt(1)-1,2,'VISp','color',color1(1,:),'fontSize',16,'Interpreter','none');
% text(qt(1)-1,4,'SSp-ul','color',color1(2,:),'fontSize',16,'Interpreter','none');
% text(qt(1)-1,6,'SSp-m','color',color1(3,:),'fontSize',16,'Interpreter','none');
xlabel('Time (s)','fontSize',16);