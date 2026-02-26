function h1ac = plotSpiralTimeSeries(data_folder,save_folder)
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
% load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
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
params.downscale = 1;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
%%
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
freq = [2,8];
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
Ur = reshape(Utransformed, size(Utransformed,1)*size(Utransformed,2), size(Utransformed,3));
% rawTrace1 = Ur*dV1;
rawTrace1 = Ur*V(:,frameTemp);
rawTrace1 = rawTrace1 -mean(rawTrace1 ,2);
tsize = size(dV1,2);
tq = 1:rate:tsize;
rawTrace1 = interp1(1:tsize,rawTrace1',tq);
rawTrace1 = rawTrace1';
rawTrace1 = reshape(rawTrace1,size(Utransformed,1),size(Utransformed,2),[]);
rawTrace = rawTrace1(:,:,1+35/rate:end-35/rate);
rawTrace = rawTrace./mimgtransformed;
%%
t1 = t(frameStart:frameEnd);

%%

color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
frames_to_plot = 18;
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
%%
scale3 = 5;
h1ac = figure('Renderer', 'painters', 'Position', [100 100 800 800]);
ax1 = subplot(1,1,1);
frame = 21;
framea = tracePhase1(:,:,frame);
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; 
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
%%
% Your phase matrix (values from -pi to pi)
phaseData = framea; % your 2D matrix
% phaseData(~BW2) = nan;
% phaseData = reshape(phaseData, size(framea,1),size(framea,2));
% Normalize to [0, 1]
normData = (phaseData + pi) / (2 * pi);
% Convert to colormap indices (e.g., 256-level)
nColors = 256;
ind = round(normData * (nColors - 1)) + 1;
ind = min(max(ind, 1), nColors); % clamp
% Get the colormap (e.g., hsv, parula, jet, etc.)
cmap = colorcet('C06','N',nColors);
% % Convert to RGB image (H x W x 3)
rgbImage = ind2rgb(ind, cmap);
% rgbimage2 = imresize(rgbImage, [1320,1140]);
figure;
image(rgbImage)
%%
save('phase_colormap.mat', 'rgbImage');
% writeNPY('rgbImage', 'phase_colormap.npy')
%%
print(h1ac,fullfile(save_folder, 'Fig1ac_example_time_series&flow_V.pdf'),...
    '-dpdf', '-bestfit', '-painters');