function h1ac = instantenous_frequency_example1(data_folder,save_folder)
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
rate = 0.1;
%%
% pixel(1,:) = [900,800]; % VISp
% pixel(2,:) = [775,650]; % RSP
% pixel(3,:) = [590,750]; % SSp-ul
% pixel(4,:) = [520,850]; % SSp-ll
% pixel(5,:) = [480,960]; % SSp-m
% pixel(6,:) = [550,960]; % SSp-n
% pixel(7,:) = [682,950]; % SSp-bfd
% pixel = round(pixel/params.downscale);
%%
center = [79,110];
th2 = 1:45:360; 
px1 = center(2);
py1 = center(1);
r = 12;
cx2 = round(r*cosd(th2)+px1);
cy2 = round(r*sind(th2)+py1);
pixel = [cy2;cx2]';
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
freq = [2,6];
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
tracePhase1a = unwrap(tracePhase1,[],3);
phase_diff = abs(diff(tracePhase1a,1,3));
ft_all = phase_diff*35./(2*pi);
%%
Ur = reshape(Utransformed, size(Utransformed,1)*size(Utransformed,2), size(Utransformed,3));
rawTrace1 = Ur*dV1;
rawTrace1 = rawTrace1 -mean(rawTrace1 ,2);
tsize = size(dV1,2);
tq = 1:rate:tsize;
rawTrace1 = interp1(1:tsize,rawTrace1',tq);
rawTrace1 = rawTrace1';
rawTrace = reshape(rawTrace1,size(Utransformed,1),size(Utransformed,2),[]);
rawTrace = rawTrace(:,:,1+35/rate:end-35/rate);
rawTrace = rawTrace./mimgtransformed;
%%
trace_filt = zeros(3,size(trace2d1,3));
trace_raw = zeros(3,size(trace2d1,3));
trace_phase = zeros(3,size(trace2d1,3));
for i = 1:7
    trace_raw(i,:) = rawTrace(pixel(i,1),pixel(i,2),:);
    trace_filt(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
    trace_phase(i,:) = tracePhase1(pixel(i,1),pixel(i,2),:);
end
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
scale3 = 5/8;
% color2 = cbrewer2('seq','YlOrRd',9);
color2 = cbrewer2('seq','YlOrRd',10);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
frames_to_plot = 18;
first_frame =22;
last_frame = first_frame+frames_to_plot;
t1a = t1(first_frame); t1b = t1(last_frame);
first_frame1 = find(qt-t1a>0, 1, 'first');
last_frame1 = find(qt-t1b>0, 1, 'first');

example_frame_to_plot = first_frame+8;
t1c = t1(example_frame_to_plot); 
frame = find(qt-t1c>0, 1, 'first');
dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
%%
color1 = cbrewer2('div','RdBu',20);
color1 = flipud(color1);
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/8;
frames = 10;
subplotn = 10;
first_frame = 23;
% first_frame = 160;

dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
th2 = 1:5:360; 
h1 = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+(k-1);
    frame_real = frameStart+first_frame+k-1;
    %%
    ax3 = subplottight(3,subplotn,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    % caxis([raw_min*0.9,raw_max*0.9]);
    caxis([-1.5,1.5]);
    %%    
    ax4 = subplottight(3,subplotn,subplotn+k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    %%
    ax4 = subplottight(3,subplotn,subplotn*2+k);
    im_phase = imagesc(ft_all(:,:,frame));
    colormap(ax4,color1);
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([0,8]);
    frame_count = frame_count+1;
end
ax3 = subplottight(3,subplotn,k);
im_raw = imagesc(dff(:,:,frame));
colormap(ax3, parula)
axis image; axis off;
set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([-1.5,1.5]);
cb1 = colorbar;
ax4 = subplottight(2,subplotn,subplotn+k);
im_phase = imagesc(ft_all(:,:,frame));
colormap(ax4,color1);
axis image; axis off;
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([0,8]);
cb2 = colorbar;
cb2.Ticks = [0:4:8];
cb2.TickLabels = num2cell([0:4:8]);

print(h1, 'instananeous frequency_example1.pdf',...
    '-dpdf', '-bestfit', '-painters');