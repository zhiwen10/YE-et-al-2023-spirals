function plotSpiralSequence_freq2(T,mouseID,freq,data_folder,save_folder)
%% create save folder
freq_name = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
freq_folder = fullfile(save_folder,freq_name);
if ~exist(freq_folder, 'dir')
    mkdir(freq_folder);
end
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% load data from an example an session
id = find(strcmp(T.MouseID,mouseID));
mn = T.MouseID{id};
tda = T.date(id);
en = T.folder(id);    
td = datestr(tda,'yyyy-mm-dd');
    
tdb = datestr(td,'yyyymmdd');
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals\svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
%% registration
fname = [mn '_' tdb '_' num2str(en)];
% load(fullfile(data_folder,'tables',[fname '_tform.mat']));               % load atlas transformation matrix tform;
load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));
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
% yscaling = 50;                                                             % for 0.5-2Hz
yscaling = 50;                                                             % for 02-8Hz
%%
% load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
BW = logical(projectedAtlas1);
BW = BW(1:params.downscale:end,1:params.downscale:end);
BW2 = BW;
%%
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
pixel = round(pixel/params.downscale);
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
%%
clear spirals_temp spirals_temp2
clear trace2d1 tracePhase1
frameStart = 58076; 
frameEnd = 5600+frameStart;
% params.downscale = 10;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
% frameEnd =  spirals_temp(end,5)+5;
%% hilbert transform
frameTemp = frameStart-350:frameEnd+350; % extra 2*35 frames before filter data 
V1 = V(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,V1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+(350)/rate:end-(350)/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+(350)/rate:end-(350)/rate); % reduce 2*35 frames after filter data 
%%
trace = zeros(7,size(trace2d1,3));
for i = 1:7
    trace(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
end
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
frames = 36;
col = 18;
first_frame =1;
dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
%%
h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+140*(k-1);
    ax3 = subplottight(4,col,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([raw_min*0.6,raw_max*0.6]);
    text_ax1 = title(['frame',num2str(frame_count)],'fontSize',8);

    ax4 = subplottight(4,col,36+k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'))
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    frame_count = frame_count+1;
end
print(h1b, fullfile(freq_folder,...
    ['Spiral_sequence_' num2str(frameStart) '.pdf']),...
    '-dpdf', '-bestfit', '-painters');
%%
rate1 = 1;
[trace2d1,traceAmp1,tracePhase1,rawTrace1] = spiralPhaseMap5_freq(Utransformed,V1,t,params,freq,rate1);
trace2d1 = trace2d1(:,:,1+350/rate1:end-350/rate1)./mimgtransformed; 
rawTrace = rawTrace1(:,:,1+350/rate1:end-350/rate1)./mimgtransformed; 
dff = trace2d1*100;
tracePhase1 = tracePhase1(:,:,1+350/rate1:end-350/rate1); % reduce 2*35 frames after filter data 

trace_filt = zeros(7,size(trace2d1,3));
trace_raw = zeros(7,size(trace2d1,3));
trace_phase = zeros(7,size(trace2d1,3));
for i = 1:7
    trace_raw(i,:) = rawTrace(pixel(i,1),pixel(i,2),:);
    trace_filt(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
    trace_phase(i,:) = tracePhase1(pixel(i,1),pixel(i,2),:);
end
%
scale3 = 5/8;
t1 = t(frameStart:frameEnd);
tq = 1:rate1:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%
color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ll','SSp_ul','SSp_m','SSp_n','SSp_bfd'};
last_frame = first_frame+frames;                                       % first_frame = 1, frames = 36; defined previously
t1a = t1(first_frame); t1b = t1(last_frame);
first_frame_to_plot1 = find(qt-t1a>0,1, 'first');
last_frame_to_plot1 = find(qt-t1b>0,1, 'first');

raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -3.14; phase_max = 3.14;

h2 = figure('Renderer', 'painters', 'Position', [100 100 700 500]);
ax1 = subplot(1,2,1);
im1= imshow(mimgtransformedRGB);
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(3:end,:),'filled');

ax2 = subplot(1,2,2);
for i = 1:7
    plot(ax2,qt,trace_raw(i,:)*yscaling+3*(i-1),'k','lineWidth',1);
    hold on; 
    color1a = trace_phase(i,:)/(2*pi)+0.5;
    color_hsv = colormap(ax2,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    color_hsv1 = interp1(num1,color_hsv,color1a);
    qt1 = qt(1:end-1);qt2 = qt(2:end);
    trace1 = trace_filt(i,1:end-1);
    trace2 = trace_filt(i,2:end);
    for j = 1:numel(qt1)
    plot(ax2,[qt1(j),qt2(j)],[trace1(j)*yscaling+3*(i-1),trace2(j)*yscaling+3*(i-1)],...
        'lineWidth',2,'color',color_hsv1(j,:));
    hold on; 
    end
end
ax2.FontSize = 8; 
xlim([t1(1),t1(end)]);

set(ax2,'YTickLabel',[]);
for i = 1:7
    text(t1(1)-0.5,(i-1)*3,nameList{i},'Color', color2(i+2,:),'Interpreter','None');
end
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(qt(first_frame_to_plot1),'--');
line_ax2 = xline(qt(last_frame_to_plot1),'--');
hold on;
plot([t1(1)+0.5,t1(1)+0.5],[0,2],'r');
print(h2, fullfile(freq_folder,...
    ['Spiral_sequence_' num2str(frameStart) '-2.pdf']),...
    '-dpdf', '-bestfit', '-painters');
close all;
