function [hs1i,hs1j] = instantenous_frequency_example3(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
kk = 12;                                                                   % ZYE_0067
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform;
load(fullfile(data_folder,'spirals','spirals_example',[fname '_mask.mat']));
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed1 = Utransformed(1:8:end,1:8:end,:);
mimgtransformed = mimgtransformed(1:8:end,1:8:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
first_frame =35;
frameStart = 89053-first_frame; frameEnd = frameStart +70;
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1,rawTrace1] = spiralPhaseMap5(Utransformed1,dV1,t,params,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
rawTrace = rawTrace1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
dff = trace2d1*100;
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
tracePhase1a = unwrap(tracePhase1,[],3);
phase_diff = abs(diff(tracePhase1a,1,3));
ft_all = phase_diff*35./(2*pi);
%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ll','SSp_ul','SSp_m','SSp_n','SSp_bfd'};
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,920]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
pixel(:,2) = 1140-pixel(:,2);
pixel = round(pixel/params.downscale);
%%
trace_filt = zeros(7,size(trace2d1,3));
trace_raw = zeros(7,size(trace2d1,3));
trace_phase = zeros(7,size(trace2d1,3));
for i = 1:7
    trace_raw(i,:) = rawTrace(pixel(i,1),pixel(i,2),:);
    trace_filt(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
    trace_phase(i,:) = tracePhase1(pixel(i,1),pixel(i,2),:);
end
scale3 = 5/8;
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
a = 45; b = 95; % height
c = 20; d = 70; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
%%
color1 = cbrewer2('div','RdBu',20);
color1 = flipud(color1);
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/8;
last_frame = first_frame+18;
rows = 1;
raw_min =-2; raw_max = 2;

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
    caxis([raw_min,raw_max]);
    % caxis([-1.5,1.5]);
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
caxis([raw_min,raw_max]);
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
%%
print(h1, 'instananeous frequency_example3.pdf',...
    '-dpdf', '-bestfit', '-painters');