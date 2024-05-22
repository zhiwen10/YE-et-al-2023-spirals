function plotFftnVideo(T,data_folder,save_folder)
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
load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));       % load atlas transformation matrix tform
load(fullfile(data_folder,'spirals\fft_roi',[fname '_roi.mat']));          % read rectangle ROIouseID{kk};%% read svd components (U,V,t) from processed data folder
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals\svd',subfolder);
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
tStart = 1681; tEnd = 1685;                                                % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); 
frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35;                                     % extra 2*35 frames before filter data 
dV1 = dV(1:50,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = ...
    spiralPhaseMap_freq(U1,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate);                      % reduce 2*35 frames after filter data 
%% 3d fft
[trace2d1_fft,traceAmp1_fft,tracePhase1_fft] = ...
    spiralPhaseMap_fftn_freq(U1,dV1,t,params,freq,rate);
trace2d1_fft = trace2d1_fft(:,:,1+35/rate:end-35/rate); 
tracePhase1_fft = tracePhase1_fft(:,:,1+35/rate:end-35/rate);              % reduce 2*35 frames after filter data 
%%
dff = trace2d1*100;
dff_fft = trace2d1_fft*100;
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
filename = fullfile(save_folder,fname);
v = VideoWriter(filename);
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
cb1 = colorbar; 
cb1.FontSize = 16;
cb1.Color = 'w';
pos1 = get(cb1,'Position');
set(cb1,'position',[pos1(1)+0.08 pos1(2)+0.1 pos1(3) pos1(4)/2]);
caxis([raw_min,raw_max]);
text_ax1 = text(0,-20,[num2str(round(qt(1)*10)/10) ' s'],...
    'fontSize',16,'color','w');
text1_ax1 = text(150,-50,'Raw trace (df/f, %)',...
    'fontSize',16,'color','w');

ax2 = subplot(2,2,2);
im_phase = imagesc(tracePhase1(:,:,1));
colormap(ax2,colorcet('C06'));
axis image; axis off;
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
cb1 = colorbar; 
cb1.FontSize = 16;
cb1.Color = 'w';
pos1 = get(cb1,'Position');
set(cb1,'position',[pos1(1)+0.08 pos1(2)+0.1 pos1(3) pos1(4)/2]);
caxis([raw_min,raw_max]);
text3_ax1 = text(150,-50,'3D-FFT trace (df/f, %)',...
    'fontSize',16,'color','w');

ax4 = subplot(2,2,4);
im_phase_fft = imagesc(tracePhase1_fft(:,:,1));
colormap(ax4,colorcet('C06'));
axis image; axis off;
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
[trace2d1a,traceAmp1a,tracePhase1a] = ...
    spiralPhaseMap_freq(U(:,:,1:50),dV1,t,params,freq,rate1);
trace2d1a = trace2d1a(:,:,1+35/rate1:end-35/rate1); 
tracePhase1a = tracePhase1a(:,:,1+35/rate1:end-35/rate1);                  % reduce 2*35 frames after filter data 
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
patch('Faces',f1,'Vertices',v1,'FaceColor','None',...
    'EdgeColor','k','lineWidth',2);
print(h1, fullfile(save_folder,'ZYE12_FFT_ROI.pdf'),...
    '-dpdf', '-bestfit', '-painters');
