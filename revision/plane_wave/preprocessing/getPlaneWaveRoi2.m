function getPlaneWaveRoi2(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
params.downscale = 8;
%%
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
kk = 5;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals\svd',subfolder);
mimg = readNPY(fullfile(session_root, 'meanImage.npy'));
% registration
load(fullfile(data_folder,'spirals\rf_tform_8x',[fname '_tform_8x.mat']));   % load atlas transformation matrix tform;
% get flow field from unregistered frame frist, then transform to
% registered, to avoid interpolation problem with circular phase 
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
mimgt = imwarp(mimg1,tform,'OutputView',imref2d(size(BW1)));
%% save boarder rectangle roi
h1 = figure('Renderer', 'painters', 'Position', [50 50 800 350]);
ax1 = subplot(1,1,1);
im1 = imagesc(mimgt);
colormap(ax1,'gray');
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
save(fullfile(save_folder,'border_roi2.mat'),'roi','bwroi');
