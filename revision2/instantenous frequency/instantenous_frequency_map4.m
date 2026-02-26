function h1b = instantenous_frequency_map4(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% load data from an example an session
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
load(fullfile(data_folder,'spirals','spirals_grouping',...
    [fname '_spirals_group_fftn.mat']));
% load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat'])); 
spirals = cell2mat(archiveCell);   
clear spiralsT 
[spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
    tform,spirals(:,1),spirals(:,2));    
spirals(:,1:2) = round(spiralsT/8); 
%%
pixSize = 3.45/1000/0.6*3;  % mm / pix for spiral radius
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
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
pixel(:,2) = 1140-pixel(:,2);
pixel = round(pixel/params.downscale);
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed1 = Utransformed(1:8:end,1:8:end,:);
mimgtransformed = mimgtransformed(1:8:end,1:8:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%% hilbert transform
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
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
figure;
for i = 1:7
    subplot(7,1,i);
    % plot(qt(1:end-1),squeeze(phase_diff(pixel(i,1),pixel(i,2),:))*35/(2*pi));
    plot(qt(1:end-1),squeeze(ft_all(pixel(i,1),pixel(i,2),:)));
    hold on;
    plot(qt,squeeze(tracePhase1(pixel(i,1),pixel(i,2),:)));
    ylim([-3.14,8]);
end
%%
tracePhaseDiff = diff(tracePhase1,1,3);
ft_all = angle(exp(i*(tracePhaseDiff)))*35./(2*pi);
figure;
for i = 1:7
    subplot(7,1,i);
    plot(qt(1:end-1),squeeze(ft_all(pixel(i,1),pixel(i,2),:)));
    hold on;
    plot(qt,squeeze(tracePhase1(pixel(i,1),pixel(i,2),:)));
    ylim([0,8]);
end
%%
figure;
for i = 1:7
    subplot(7,1,i);
    phase1 = squeeze(tracePhase1(pixel(i,1),pixel(i,2),:));
    freq1 = squeeze(ft_all(pixel(i,1),pixel(i,2),:));
    scatter(phase1(1:end-1),freq1);
    ylim([0,6]);
end
%%
trace = zeros(7,size(trace2d1,3));
for k = 1:7
    trace(k,:) = trace2d1(pixel(k,1),pixel(k,2),:);
end

%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
scale3 = 5/8;
%%
% a = 50; b = 110; % height
% c = 75; d = 135; % width
a = 63; b = 100; % height
c = 95; d = 132; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];

lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/8;
frames = 20;
subplotn = 10;
first_frame = 23;
% first_frame = 160;

dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
th2 = 1:5:360; 

h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+(k-1);
    frame_real = frameStart+first_frame+k-1;
    %%    
    ax4 = subplottight(4,subplotn,k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    %%
    ax4 = subplottight(4,subplotn,subplotn*2+k);
    im_phase = imagesc(ft_all(:,:,frame));
    colormap(ax4,parula);
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([0,8]);

    frame_count = frame_count+1;
end
%%
print(h1b, fullfile(save_folder,'Fig1b_example_spiral_sequence3.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%%
%%
rate1 = 1;
[trace2d1,traceAmp1,tracePhase1,rawTrace1] = spiralPhaseMap5(Utransformed1,dV1,t,params,rate1);
trace2d1 = trace2d1(:,:,1+35/rate1:end-35/rate1)./mimgtransformed; 
rawTrace = rawTrace1(:,:,1+35/rate1:end-35/rate1)./mimgtransformed; 
dff = trace2d1*100;
tracePhase1 = tracePhase1(:,:,1+35/rate1:end-35/rate1); % reduce 2*35 frames after filter data 

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
rate1 = 1;
t1 = t(frameStart:frameEnd);
tq = 1:rate1:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%
color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ll','SSp_ul','SSp_m','SSp_n','SSp_bfd'};
frames_to_plot = 10;
last_frame = first_frame+frames_to_plot;
t1a = t1(first_frame); t1b = t1(last_frame);
first_frame_to_plot1 = find(qt-t1a>0,1, 'first');
last_frame_to_plot1 = find(qt-t1b>0,1, 'first');

raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;

hs1i = figure('Renderer', 'painters', 'Position', [100 100 700 500]);
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
    plot(ax2,qt,trace_raw(i,:)*100+3*(i-1),'k','lineWidth',1);
    hold on; 
    color1a = trace_phase(i,:)/(2*pi)+0.5;
    color_hsv = colormap(ax2,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    color_hsv1 = interp1(num1,color_hsv,color1a);
    qt1 = qt(1:end-1);qt2 = qt(2:end);
    trace1 = trace_filt(i,1:end-1);
    trace2 = trace_filt(i,2:end);
    for j = 1:numel(qt1)
    plot(ax2,[qt1(j),qt2(j)],[trace1(j)*100+3*(i-1),trace2(j)*100+3*(i-1)],...
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
