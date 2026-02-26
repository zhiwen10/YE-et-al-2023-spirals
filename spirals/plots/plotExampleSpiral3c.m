function [hs1i,hs1j] = plotExampleSpiral3c(T,data_folder,save_folder)
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
load(fullfile(data_folder,'spirals','spirals_grouping',...
    [fname '_spirals_group_fftn.mat']));
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat'])); 
spirals = cell2mat(archiveCell);   
clear spiralsT 
[spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
    tform,spirals(:,1),spirals(:,2));    
spirals(:,1:2) = round(spiralsT/8); 
%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ll','SSp_ul','SSp_m','SSp_n','SSp_bfd'};
% pixel(1,:) = [900,800]; % VISp
% pixel(2,:) = [775,650]; % RSP
% pixel(3,:) = [590,750]; % SSp-ul
% pixel(4,:) = [520,850]; % SSp-ll
% pixel(5,:) = [480,920]; % SSp-m
% pixel(6,:) = [550,960]; % SSp-n
% pixel(7,:) = [682,950]; % SSp-bfd
% pixel(:,2) = 1140-pixel(:,2);
% pixel = round(pixel/params.downscale);

% center = [77,36];
% r = 15;
% th2 = linspace(0,360,8)+300;

center = [74,34];
r = 13;
th2 = linspace(0,360,8)+300;
th2 = th2(1:end-1);
px1 = center(2);
py1 = center(1);

cx2 = round(r*cosd(th2)+px1);
cy2 = round(r*sind(th2)+py1);
pixel = [cy2;cx2]';
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
lineColor = 'w'; lineColor1 = 'w';
hemi = [];
last_frame = first_frame+18;
rows = 1;
pixSize = 3.45/1000/0.6*3;  % mm / pix for spiral radius

raw_min =-2; raw_max = 2;
% raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
hs1j = figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
frame_count = 1;
subplotn = 10;
th2 = 1:5:360; 
for k = 1:subplotn
    frame  =first_frame+(k-1);
    frame_real = frameStart+first_frame+k-1;
    spiral_temp = spirals(spirals(:,5) ==frame_real,:);
    ax3 = subplottight(4,subplotn,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([raw_min,raw_max]);
    text_ax1 = title(['frame',num2str(frame_count)],'fontSize',8);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    ax4 = subplottight(4,subplotn,subplotn+k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'))
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    % text_ax2 = text(0,0,[num2str(round(qt(frame)*10)/10) ' s']);
    text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    ax5(k) = subplottight(4,subplotn,subplotn*2+k);
    ax5(k).Position(1) = ax5(k).Position(1)+0.005;
    ax5(k).Position(2) = ax5(k).Position(2);
    ax5(k).Position(3) = ax5(k).Position(3)-0.015;
    ax5(k).Position(4) = ax5(k).Position(4)-0.03;
    im_phase = imagesc(dff(:,:,frame));
    colormap(ax5(k),parula)
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
%     hold on;
%     scatter(spiral_temp(:,1),spiral_temp(:,2),8,'k','filled');
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3)*pixSize/(0.01*8);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',1);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end       
    caxis([-1,1]);
    % text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    ylim([a,b]);
    xlim([c,d]);
    
    ax6(k) = subplottight(4,subplotn,subplotn*3+k);
    ax6(k).Position(1) = ax6(k).Position(1)+0.005;
    ax6(k).Position(2) = ax6(k).Position(2);
    ax6(k).Position(3) = ax6(k).Position(3)-0.015;
    ax6(k).Position(4) = ax6(k).Position(4)-0.03;
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax6(k),colorcet('C06'))
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3)*pixSize/(0.01*8);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',1);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    caxis([phase_min,phase_max]);
    % text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    ylim([a,b]);
    xlim([c,d]);
    
    frame_count = frame_count+1;
end
frameS = frameStart+first_frame;
%%
print(hs1j,fullfile(save_folder,...
    ['FigS1j_' fname '_' num2str(frameS)]),...
    '-dpdf', '-bestfit', '-painters');
%%
rate1 = 0.1;
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
rate1 = 0.1;
t1 = t(frameStart:frameEnd);
tq = 1:rate1:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%
color2 = cbrewer2('seq','YlOrRd',9);
color2 = flipud(color2(3:end,:));

% nameList = {'VISp','RSP','SSp_ll','SSp_ul','SSp_m','SSp_n','SSp_bfd'};
frames_to_plot = 9;
last_frame = first_frame+frames_to_plot;
t1a = t1(first_frame); t1b = t1(last_frame);
first_frame_to_plot1 = find(qt-t1a>0,1, 'first');
last_frame_to_plot1 = find(qt-t1b>0,1, 'first');
frame0 = first_frame_to_plot1-30;

raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;

hs1i = figure('Renderer', 'painters', 'Position', [100 100 700 500]);
ax1 = subplot(1,4,1);
im1= imshow(mimgtransformedRGB);
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(1:end,:),'filled');
hold on;
scatter(center(2),center(1),16,'k','filled');
hold on;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);

yscale = 100;
ax2 = subplot(1,4,[2,3]);
for i = 1:7
    plot(ax2,qt,trace_raw(i,:)*yscale-3*(i-1),'k','lineWidth',1);
    hold on; 
    color1a = trace_phase(i,:)/(2*pi)+0.5;
    color_hsv = colormap(ax2,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    color_hsv1 = interp1(num1,color_hsv,color1a);
    qt1 = qt(1:end-1);qt2 = qt(2:end);
    trace1 = trace_filt(i,1:end-1);
    trace2 = trace_filt(i,2:end);
    for j = 1:numel(qt1)
    plot(ax2,[qt1(j),qt2(j)],[trace1(j)*yscale-3*(i-1),trace2(j)*yscale-3*(i-1)],...
        'lineWidth',2,'color',color_hsv1(j,:));
    hold on; 
    end   
    %
    trace_filt2 = trace_filt(i,frame0:last_frame_to_plot1+10);
    [pks,locs] = findpeaks(trace_filt2,'MinPeakProminence',0.001);
    scatter(qt1(frame0+locs),trace_filt2(locs)*yscale-3*(i-1)+0.2,12,'k','filled');    
    hold on; 
end
ax2.FontSize = 8; 
xlim([t1(1),t1(end)]);
set(ax2,'YTickLabel',[]);
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(qt(first_frame_to_plot1),'--');
line_ax2 = xline(qt(last_frame_to_plot1),'--');
hold on;
plot([t1(1)+0.5,t1(1)+0.5],[0,2],'r');

qta = qt-qt(first_frame_to_plot1);
ax2 = subplot(1,4,4);
for i = 1:7
    plot(ax2,qta,trace_raw(i,:)*yscale-3*(i-1),'k','lineWidth',1);
    hold on; 
    color1a = trace_phase(i,:)/(2*pi)+0.5;
    color_hsv = colormap(ax2,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    color_hsv1 = interp1(num1,color_hsv,color1a);
    qt1 = qta(1:end-1);qt2 = qta(2:end);
    trace1 = trace_filt(i,1:end-1);
    trace2 = trace_filt(i,2:end);
    for j = 1:numel(qt1)
    plot(ax2,[qt1(j),qt2(j)],[trace1(j)*yscale-3*(i-1),trace2(j)*yscale-3*(i-1)],...
        'lineWidth',2,'color',color_hsv1(j,:));
    hold on; 
    end
    
    trace_filt2 = trace_filt(i,frame0:last_frame_to_plot1+10);
    [pks,locs] = findpeaks(trace_filt2,'MinPeakProminence',0.001);
    scatter(qt1(frame0+locs),trace_filt2(locs)*yscale-3*(i-1)+0.2,12,'k','filled');    
    hold on; 
end
ax2.FontSize = 8; 
xlim([qta(frame0),qta(last_frame_to_plot1+10)]);
set(ax2,'YTickLabel',[]);
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(qta(first_frame_to_plot1),'--');
line_ax2 = xline(qta(last_frame_to_plot1),'--');
hold on;
plot([t1(1)+0.5,t1(1)+0.5],[0,2],'r');
%%
print(hs1i, fullfile(save_folder,...
    ['FigS1i_' fname '_' num2str(frameS) '_cetc6_4.pdf']),...
    '-dpdf', '-bestfit', '-painters');