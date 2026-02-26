function [h1b,r_var] = plotSpiralSequences2(U_predict,dV_predict,U_raw,dV_raw,t,centers,radius,frames)
%%
data_folder1 = 'E:\spiral_data_share\data';
load(fullfile(data_folder1,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder1,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder1);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%%
freq = [2,8];
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
frameTemp = frames(1)-35:frames(end)+35; % extra 2*35 frames before filter data 
dV1 = dV_predict(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1_pred] = spiralPhaseMap_freq(U_predict,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1_pred = tracePhase1_pred(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
t1 = t(frames);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
scale3 = 5/8;
%%
BW2 = BW(1:8:end,1:8:end);
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
th2 = 1:5:360; 
framesa = numel(frames);
if framesa>9
    framesa = 9;
end
subplotn = 10;
first_frame =1;
dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
frame_count = 1;
for k = 1:framesa
    frame  =first_frame+(k-1);
    frame_real = frames(1)+first_frame+k-1;
    %%
    ax3 = subplottight(5,subplotn,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([raw_min,raw_max]);
    text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    
    for kk = 1:size(centers,1)
        px1 = centers(kk,1);
        py1 = centers(kk,2);
        r = radius(kk,1);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
    
    ax4 = subplottight(5,subplotn,subplotn+k);
    im_phase = imagesc(tracePhase1_pred(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);   
    
    for kk = 1:size(centers,1)
        px1 = centers(kk,1);
        py1 = centers(kk,2);
        r = radius(kk,1);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
      
    frame_count = frame_count+1;
end
ax4 = subplottight(5,subplotn,10);
im_raw = imagesc(dff(:,:,frame));
colormap(ax3, parula)
axis image; axis off;
set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([raw_min,raw_max]);
text_ax1 = title(['frame',num2str(frame_count)],'fontSize',8);
colorbar
%%
frameTemp = frames(1)-35:frames(end)+35; % extra 2*35 frames before filter data 
dV1 = dV_raw(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1_raw] = spiralPhaseMap_freq(U_raw,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate); 
tracePhase1_raw = tracePhase1_raw(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
dff = trace2d1*100;
% raw_min = min(dff(:)); raw_max = max(dff(:));
% phase_min = -pi; phase_max = pi;
frame_count = 1;
for k = 1:framesa
    frame  =first_frame+(k-1);
    frame_real = frames(1)+first_frame+k-1;
    %%
    ax3 = subplottight(5,subplotn,k+2*subplotn);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([raw_min,raw_max]);
    
    for kk = 1:size(centers,1)
        px1 = centers(kk,1);
        py1 = centers(kk,2);
        r = radius(kk,1);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
    
    ax4 = subplottight(5,subplotn,3*subplotn+k);
    im_phase = imagesc(tracePhase1_raw(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    
    for kk = 1:size(centers,1)
        px1 = centers(kk,1);
        py1 = centers(kk,2);
        r = radius(kk,1);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
      
    frame_count = frame_count+1;
end
ax4 = subplottight(5,subplotn,30);
im_raw = imagesc(dff(:,:,frame));
colormap(ax3, parula)
axis image; axis off;
set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([raw_min,raw_max]);
text_ax1 = title(['frame',num2str(frame_count)],'fontSize',8);
colorbar

%%
phase_diff = angdiff(tracePhase1_pred,tracePhase1_raw);
r_var = circ_r(phase_diff, [], [], 3);
for k = 1:framesa
    frame  =first_frame+(k-1);
    ax4 = subplottight(5,subplotn,4*subplotn+k);
    im_phase = imagesc(phase_diff(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    
    for kk = 1:size(centers,1)
        px1 = centers(kk,1);
        py1 = centers(kk,2);
        r = radius(kk,1);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
      
    frame_count = frame_count+1;
end
