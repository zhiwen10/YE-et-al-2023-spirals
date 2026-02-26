function h1ac = plotSpiralTimeSeries3(data_folder,save_folder)
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
load(fullfile(data_folder,'tables','mask_ZYE12_3.mat'));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'tables',[fname '_tform_2.mat']));                 % load atlas transformation matrix tform;
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
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
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
freq = [2,8];
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
% tStart = 1660; tEnd = 1663; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
Ur = reshape(Utransformed, size(Utransformed,1)*size(Utransformed,2), size(Utransformed,3));
rawTrace1 = Ur*dV(:,frameTemp);
rawTrace1 = rawTrace1 -mean(rawTrace1 ,2);
tsize = size(dV1,2);
tq = 1:rate:tsize;
rawTrace1 = interp1(1:tsize,rawTrace1',tq);
rawTrace1 = rawTrace1';
rawTrace1 = reshape(rawTrace1,size(Utransformed,1),size(Utransformed,2),[]);
rawTrace = rawTrace1(:,:,1+35/rate:end-35/rate);
rawTrace = rawTrace./mimgtransformed;
%%
% params1
% center = [78,112];
% r = 15;
% params2
% center = [77,112];
% r = 15;
% params3
% center = [77,113];
% r = 17;
% params4
center = [75,111];
r = 15;
th2 = linspace(0,360,8)+294;
th2 = th2(1:end-1);
px1 = center(2);
py1 = center(1);

cx2 = round(r*cosd(th2)+px1);
cy2 = round(r*sind(th2)+py1);
pixel = [cy2;cx2]';
%%
trace_filt = zeros(3,size(trace2d1,3));
trace_raw = zeros(3,size(trace2d1,3));
trace_phase = zeros(3,size(trace2d1,3));
for i = 1:size(pixel,1)
    trace_raw(i,:) = rawTrace(pixel(i,1),pixel(i,2),:);
    trace_filt(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
    trace_phase(i,:) = tracePhase1(pixel(i,1),pixel(i,2),:);
end
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%
scale3 = 5/8;
% color2 = cbrewer2('seq','YlOrRd',9);
color2 = cbrewer2('seq','YlOrRd',size(pixel,1)+2);
% nameList = {'SSp_tr','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
nameList = {'SSp_ul','SSp_tr','SSp_bfd','SSp_n','SSp_m','SSp_ll'};
order1 = [1,2,3,3,3,4,6,6];
frames_to_plot = 9;
first_frame = 22;
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
%
h1ac = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
ax1 = subplot(2,2,1);
im1= imshow(mimgtransformedRGB);
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);

cb1 = colorbar(ax1);
pos1 = cb1.Position;
cb1.Position = [pos1(1)+0.1, pos1(2),pos1(3),pos1(4)/2];
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(3:end,:),'filled');
hold on;
scatter(center(2),center(1),16,'k','filled');

ax4 = subplot(2,2,2);
framea = tracePhase1(:,:,frame);
frameb = tracePhase1(:,:,frame+10);
frame_ab(1,:,:) = framea;
frame_ab(2,:,:) = frameb;
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(frame_ab,useGPU);
vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
skip = 3;
zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW2) = nan; vyRaw2(~BW2) = nan;

hold off;
im_phase = imagesc(framea);
colormap(ax4,colorcet('C06'));
axis image; 
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
scatter(center(2),center(1),16,'k','filled');
v = [95 60; 125 60;125,100;95 100];
% a = 63; b = 100; % height
% c = 95; d = 132; % width
% % rect_roi1 = [a:b]; % height
% % rect_roi2 = [c:d]; % width
% v = [c a;d a; d b;c b];
f = [1,2,3,4];
patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',2);

cb2 = colorbar(ax4);
pos2 = cb2.Position;
cb2.Position = [pos2(1)+0.1, pos2(2),pos2(3),pos2(4)/2];

ax2 = subplot(2,2,3);
yscale = 150;
for i = 1:size(pixel,1)
    plot(ax2,qt,trace_raw(i,:)*yscale+3*(i-1),'k','lineWidth',1);
    hold on; 
    color1 = trace_phase(i,:)/(2*pi)+0.5; % from [-0.5,0.5] to [0,1]
    color_hsv = colormap(ax2,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    color_hsv1 = interp1(num1,color_hsv,color1);
    qt1 = qt(1:end-1);qt2 = qt(2:end);
    trace1 = trace_filt(i,1:end-1);
    trace2 = trace_filt(i,2:end);
    for j = 1:numel(qt1)
    plot(ax2,[qt1(j),qt2(j)],[trace1(j)*yscale+3*(i-1),trace2(j)*yscale+3*(i-1)],...
        'lineWidth',2,'color',color_hsv1(j,:));
    hold on; 
    end
end
ax2.FontSize = 8; 
xlim([t1(1),t1(end)]);
set(ax2,'YTickLabel',[]);
% for i = 1:size(pixel,1)
%     text(t1(1)-0.8,(i-1)*3,nameList{order1(i)},'Color', color2(i+2,:),'Interpreter','None');
% end
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(qt(first_frame1),'--');
line_ax2 = xline(qt(last_frame1),'--');
hold on;
plot([t1(1)+0.5,t1(1)+0.5],[0,1.5],'r');
hold on;
for i = 1:size(pixel,1)
    scatter(qt1(1),3*(i-1),16,color2(2+i,:),'filled');
    hold on;
end
%
framea_zoom = framea(60:100,95:125);
ax5 = subplot(2,2,4);
hold off;
im_phase = imagesc(framea_zoom);
colormap(ax5,colorcet('C06'));
axis image; 
axis off;
hold on;
imH1Raw4 = quiver(vxRaw2(60:100,95:125),vyRaw2(60:100,95:125),'k','lineWidth',1,'autoScale','off');
cb5 = colorbar(ax5);
pos5 = cb5.Position;
cb5.Position = [pos5(1)+0.1, pos5(2),pos5(3),pos5(4)/2];
%%
print(h1ac,fullfile(save_folder, 'Fig1ac_example_time_series&flow_3.pdf'),...
    '-dpdf', '-bestfit', '-painters');