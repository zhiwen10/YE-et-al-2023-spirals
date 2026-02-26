function h1b = plotSpiralSequenceMask(data_folder,save_folder)
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
load(fullfile(data_folder,'tables',[fname '_tform.mat']));               % load atlas transformation matrix tform;
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
pixel = round(pixel/params.downscale);
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%% hilbert transform
freq = [2,8];
tStart = 1681; tEnd = 1686; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
traceAmp1 = traceAmp1(:,:,1+35/rate:end-35/rate)./mimgtransformed;
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
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
scale3 = 5/8;
%%
a = 63; b = 100; % height
c = 95; d = 132; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
BW3 = zeros(size(BW2));
BW3(a:b,c:d) = 1;
BW3 = logical(BW3);
%%
traceAmp2 = zeros(size(traceAmp1));
for k = 1:size(traceAmp1,3)
    temp = squeeze(traceAmp1(:,:,k));
    temp(not(BW3)) = 0;
    traceAmp2(:,:,k) = temp;
end
traceAmp2 = traceAmp2./max(traceAmp2(:));
amp_max = max(traceAmp1(:));
%%
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
frames = 9;
subplotn = 10;
first_frame =23;

dff = trace2d1*100;
traceAmp1 = traceAmp1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
th2 = 1:5:360; 
%%
h1 = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+(k-1);
    frame_real = frameStart+first_frame+k-1;
    spiral_temp = spirals(spirals(:,5) ==frame_real,:);
    %%
    ax1 = subplottight(6,subplotn,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax1, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([-1.5,1.5]);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);    
    
    ax2 = subplottight(6,subplotn,subplotn+k);
    im_raw = imagesc(tracePhase1(:,:,frame));
    colormap(ax2,colorcet('C06'))
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
    
    ax3 = subplottight(6,subplotn,subplotn*2+k);
    im_raw = imagesc(traceAmp1(:,:,frame));
    colormap(ax3, gray)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([0,1.5]);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
    
    ax4(k) = subplottight(6,subplotn,subplotn*3+k);
    ax4(k).Position(1) = ax4(k).Position(1)+0.005;
    ax4(k).Position(2) = ax4(k).Position(2);
    ax4(k).Position(3) = ax4(k).Position(3)-0.015;
    ax4(k).Position(4) = ax4(k).Position(4)-0.03;
    im_phase = imagesc(traceAmp1(:,:,frame));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    colormap(ax4(k),gray)
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3)*pixSize/(0.01*8);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                              % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    caxis([0,1.5]);
    ylim([a,b]);
    xlim([c,d]);
    
    ax5(k) = subplottight(6,subplotn,subplotn*4+k);
    ax5(k).Position(1) = ax5(k).Position(1)+0.005;
    ax5(k).Position(2) = ax5(k).Position(2);
    ax5(k).Position(3) = ax5(k).Position(3)-0.015;
    ax5(k).Position(4) = ax5(k).Position(4)-0.03;
    im_phase = imagesc(tracePhase1(:,:,frame));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    colormap(ax5(k),colorcet('C06'))
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3)*pixSize/(0.01*8);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                              % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    caxis([phase_min,phase_max]);
    ylim([a,b]);
    xlim([c,d]);
        
    ax6(k) = subplottight(6,subplotn,subplotn*5+k);
    ax6(k).Position(1) = ax6(k).Position(1)+0.005;
    ax6(k).Position(2) = ax6(k).Position(2);
    ax6(k).Position(3) = ax6(k).Position(3)-0.015;
    ax6(k).Position(4) = ax6(k).Position(4)-0.03;
    im_phase = imagesc(tracePhase1(:,:,frame));
    axis image; axis off;
    mask1 = squeeze(traceAmp2(:,:,frame));
    set(im_phase, 'AlphaData', mask1, 'AlphaDataMapping', 'scaled');
    hold on;
    colormap(ax6(k),colorcet('C06'))
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3)*pixSize/(0.01*8);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                              % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    caxis([phase_min,phase_max]);
    ylim([a,b]);
    xlim([c,d]);
    
    frame_count = frame_count+1;
end
ax1 = subplottight(6,subplotn,10);
im_raw = imagesc(dff(:,:,frame));
colormap(ax1, parula)
axis image; axis off;
set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
caxis([-1.5,1.5]);
hold on;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2); 
colorbar
%%
ax4(k) = subplottight(6,subplotn,subplotn*3+10);
ax4(k).Position(1) = ax4(k).Position(1)+0.005;
ax4(k).Position(2) = ax4(k).Position(2);
ax4(k).Position(3) = ax4(k).Position(3)-0.015;
ax4(k).Position(4) = ax4(k).Position(4)-0.03;
im_phase = imagesc(traceAmp1(:,:,frame));
axis image; axis off;
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
colormap(ax4(k),gray)
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
if not(isempty(spiral_temp))
    for kk = 1:size(spiral_temp,1)
        px1 = spiral_temp(kk,1);
        py1 = spiral_temp(kk,2);
        r = spiral_temp(kk,3)*pixSize/(0.01*8);
        cx2 = round(r*cosd(th2)+px1);
        cy2 = round(r*sind(th2)+py1);
        hold on;
        if spiral_temp(kk,4) == 1                                              % counterclockwise, then color white
            color1 = 'w';
        else
            color1 = 'k';                                                      % clockwise, then color black
        end
        plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
        hold on;
        scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
    end
end
caxis([0,1.5]);
ylim([a,b]);
xlim([c,d]);
colorbar;
%%    
print(h1, 'example_spiral_sequence_with_amp_mask3.pdf',...
    '-dpdf', '-bestfit', '-painters');
