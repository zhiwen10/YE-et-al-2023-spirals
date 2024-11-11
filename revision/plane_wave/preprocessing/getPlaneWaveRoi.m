function getPlaneWaveRoi(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
params.downscale = 8;
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
%%
load(fullfile(data_folder,'revision','plane_wave','example_flow_field.mat'));
%% save boarder rectangle roi
h1 = figure('Renderer', 'painters', 'Position', [50 50 800 350]);
frame = 699;
ax1 = subplot(1,1,1);
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax1,colorcet('C06'));
axis image; axis off;
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,params.downscale);
hold on;
hemi = []; scale3 = 5/8;
lineColor = 'k';
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
roi = drawpolygon(ax1);
%%
bwroi = createMask(roi);
save(fullfile(save_folder,'VISp_roi.mat'),'roi','bwroi');
