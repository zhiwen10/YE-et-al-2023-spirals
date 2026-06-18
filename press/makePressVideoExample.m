function makePressVideoExample(data_folder, save_folder)
%% atlas (use same sources as original video for exact registration match)
load('C:\Users\Steinmetz lab\OneDrive - UW\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','isocortex_horizontal_projection_outline.mat'));
[maskPath,st] = get_cortex_atlas_path(data_folder);

%% session data
mn  = 'ZYE_0012';
td  = '2020-10-16';
tdb = datestr(td,'yyyymmdd');
en  = 5;
subfolder    = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);
dV = [zeros(size(V,1),1) diff(V,[],2)];
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));

%% registration (use original tform, not tform_2)
regist_folder = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration';
fname = [mn '_' tdb '_' num2str(en) '_tform'];
load(fullfile(regist_folder, fname));

%% parameters
params.downscale = 4;
params.lowpass   = 0;
params.gsmooth   = 0;
rate      = 0.1;
% tStart    = 1680.92; tEnd = 1681.92;
tStart    = 1681.44; tEnd = 1681.88;
outScale  = 1.5;
frameRate = 30;
lineW     = 3.5;
textColor = 'k';
scale3    = 5 / params.downscale;
lineColor = 'k';
hemi      = [];
ffmpegExe = 'ffmpeg';
gifStride = 2;
gifScale  = 0.5;

%% warp to atlas and downscale
Utransformed    = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed    = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);

%% phase maps
frameStart = find(t>tStart,1,'first');
frameEnd   = find(t>tEnd,1,'first');
frameTemp  = frameStart-35:frameEnd+35;
dV1 = dV(:,frameTemp);
[~,~,tracePhase1] = spiralPhaseMap4(Utransformed,dV1,t,params,rate);
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate);

%% time vectors
t1  = t(frameStart:frameEnd);
tq  = 1:rate:numel(t1);
qt  = interp1(1:numel(t1),t1,tq);
qt1 = qt - qt(1);
nFrames = numel(qt1);

%% spatial scale bar
imgH       = size(tracePhase1,1);
imgW       = size(tracePhase1,2);
um_per_pix = 10 * params.downscale;
barLen_px  = 2000 / um_per_pix;

%% figure (axes fills frame, no padding)
ha = figure('Renderer','painters','Color','w','Units','pixels', ...
            'Position',[100 100 2*round(imgW*outScale/2) 2*round(imgH*outScale/2)]);
ha.MenuBar = 'none'; ha.ToolBar = 'none'; ha.Resize = 'off';
ax3 = axes('Parent',ha,'Position',[0 0 1 1]);

im_phase = imagesc(ax3, tracePhase1(:,:,1));
colormap(ax3, colorcet('C06'));
caxis(ax3, [-pi pi]);
axis(ax3,'image'); axis(ax3,'off');
set(ax3,'YDir','reverse');
hold(ax3,'on');

[area1] = plotOutline2(maskPath(1:11), st, atlas1, hemi, scale3, lineColor);
set(findobj(ax3,'Type','line'),'LineWidth',lineW);
BW3 = logical(imresize(BW2,[imgH,imgW]));
BW4 = imresize(area1,[imgH,imgW]);
BW3 = BW3 & BW4;
set(im_phase,'AlphaData',BW3,'AlphaDataMapping','scaled');

mgn = round(0.05*imgW);
xR  = imgW - mgn;
xL  = xR - barLen_px;
yB  = imgH - 45;
line(ax3,[xL xR],[yB yB],'Color',textColor,'LineWidth',3);
text(ax3,(xL+xR)/2,yB+5,'2 mm','Color',textColor,'FontSize',26,...
     'HorizontalAlignment','center','VerticalAlignment','top');
text_ts = text(ax3,0.03,0.95,sprintf('%.1f s',qt1(1)),...
    'Units','normalized','FontSize',26,'FontWeight','bold','Color',textColor);

%% write mp4
video_name = fullfile(save_folder,[subfolder '_' num2str(tStart) '-' num2str(tEnd)]);
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = frameRate;
v.Quality   = 100;
open(v);
for i = 1:nFrames
    set(im_phase,'CData',tracePhase1(:,:,i));
    text_ts.String = sprintf('%.1f s',qt1(i));
    writeVideo(v,getframe(ha));
end
close(v);

%% rescale mp4 to match gif output dimensions (imgW x imgH x gifScale)
mp4Raw  = [video_name '_raw.mp4'];
mp4File = [video_name '.mp4'];
movefile(mp4File, mp4Raw);
tgtW = 2*floor(imgW*gifScale);
tgtH = 2*floor(imgH*gifScale);
cmdScale = sprintf('%s -y -i "%s" -vf "scale=%d:%d:flags=lanczos" -c:v libx264 -crf 18 -pix_fmt yuv420p "%s"', ...
    ffmpegExe, mp4Raw, tgtW, tgtH, mp4File);
[sr,or] = system(cmdScale);
if sr ~= 0, warning('ffmpeg rescale failed:\n%s',or); end
delete(mp4Raw);

%% gif via ffmpeg (two-pass palette)
mp4File = [video_name '.mp4'];
gifFile = [video_name '.gif'];
palFile = [video_name '_palette.png'];
gifFps  = frameRate / gifStride;
vf = sprintf('fps=%g,scale=trunc(iw*%g/2)*2:-2:flags=lanczos',gifFps,gifScale);
cmd1 = sprintf('%s -y -i "%s" -vf "%s,palettegen=stats_mode=diff" "%s"',...
               ffmpegExe,mp4File,vf,palFile);
cmd2 = sprintf(['%s -y -i "%s" -i "%s" -lavfi '...
                '"%s[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3" "%s"'],...
               ffmpegExe,mp4File,palFile,vf,gifFile);
[s1,o1] = system(cmd1);
if s1 ~= 0
    warning('ffmpeg palettegen failed (is ffmpeg on PATH?):\n%s',o1);
else
    [s2,o2] = system(cmd2);
    if s2 ~= 0, warning('ffmpeg paletteuse failed:\n%s',o2); end
    if exist(palFile,'file'), delete(palFile); end
end
