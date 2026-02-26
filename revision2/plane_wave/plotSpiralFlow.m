function h1b = plotSpiralFlow(data_folder,save_folder)
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
frames = 9;
subplotn = 10;
first_frame =23;

dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;
th2 = 1:5:360; 

h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+(k-1);
    frame_real = frameStart+first_frame+k-1;
    spiral_temp = spirals(spirals(:,5) ==frame_real,:);
    %%
    ax3 = subplottight(4,subplotn,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([-1.5,1.5]);
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
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end 
    
    ax4 = subplottight(4,subplotn,subplotn+k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    hold on;
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
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    caxis([phase_min,phase_max]);
    %%
    ax5 = subplottight(4,subplotn,subplotn*2+k);
    framea = tracePhase1(:,:,frame);
    frameb = tracePhase1(:,:,frame+1);
    frame_ab(1,:,:) = framea;
    frame_ab(2,:,:) = frameb;
    useGPU = 0;
    [vxRaw,vyRaw] = HS_flowfield(frame_ab,useGPU);
    vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
    vxRaw2 = nan(size(vxRaw));
    vyRaw2 = nan(size(vyRaw));
    skip = 5;
    zoom_scale = 2;
    vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW2) = nan; vyRaw2(~BW2) = nan;
    
    im_phase = imagesc(framea);
    colormap(ax5,colorcet('C06'));
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    hold on;
    imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    %%
    vxy = complex(vxRaw2,vyRaw2);
    vxy_length = abs(vxy);
    vxy_length_mean = mean(vxy_length(:),'omitnan');
    %%
    spiral_temp1 = spiral_temp(spiral_temp(:,4) == 1,:);
    px1 = spiral_temp1(1,1);
    py1 = spiral_temp1(1,2);
    r = spiral_temp1(1,3)*pixSize/(0.01*8);
    [xx,yy] = meshgrid(1:size(mimgtransformed,2),1:size(mimgtransformed,1));
    pos = [xx(:),yy(:)];
    pos1 = pos-[px1,py1];
    distance1 = vecnorm(pos1,2,2);
    index1 = (distance1<r);
    map1 = zeros(size(mimgtransformed));
    map1(index1) = 1;
    vxy_spiral = vxy(logical(map1));
    vxy_length_spiral = abs(vxy_spiral);
    vxy_length_mean_spiral = mean(vxy_length_spiral(:),'omitnan');
    %%
    ax6 = subplottight(4,subplotn,subplotn*3+k); 
    im_phase = imagesc(map1);
    % colormap(ax5,colorcet('C06'));
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    %%   
    frame_count = frame_count+1;
end
%%
print(h1b, fullfile(save_folder,'Fig1b_example_spiral_sequence3.pdf'),...
    '-dpdf', '-bestfit', '-painters');
