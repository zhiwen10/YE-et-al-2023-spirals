function hs8p = plotBorderPlaneWave(data_folder,save_folder)
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
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
kk = 5;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root); 
dV = [zeros(size(V,1),1) diff(V,[],2)];
% registration
load(fullfile(data_folder,'spirals','rf_tform_8x',[fname '_tform_8x.mat']));   % load atlas transformation matrix tform;
% get flow field from unregistered frame frist, then transform to
% registered, to avoid interpolation problem with circular phase 
U1 = U(1:params.downscale:end,1:params.downscale:end,:);
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
U1 = U1./mimg1;
epoch = 990*35:992*35;
frame = 11;
[~,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1(:,:,1:50),dV(1:50,epoch),t,params,freq,rate);
% clear U U1 V dV
tracePhase1 = permute(tracePhase1,[3,1,2]);
useGPU = 1;
[vxRaw,vyRaw] = HS_flowfield(tracePhase1,useGPU);
% now let's transform flow field to registered atlas
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
vxRaw1 = permute(vxRaw,[2,3,1]);
vyRaw1 = permute(vyRaw,[2,3,1]);
vxRawt = imwarp(vxRaw1,tform,'OutputView',imref2d(size(BW1)));
vyRawt = imwarp(vyRaw1,tform,'OutputView',imref2d(size(BW1)));
% get phase map from artlas transformed U space
Ut = imwarp(U1(:,:,1:50),tform,'OutputView',imref2d(size(BW1)));
[~,traceAmp1t,tracePhase1t] = spiralPhaseMap_freq(Ut,dV(1:50,epoch),t,params,freq,rate);
tracePhase1t = permute(tracePhase1t,[3,1,2]);
%%
hs8p = figure('Renderer', 'painters', 'Position', [50 50 900 700]);
frame = 11;
lineColor = 'k';
hemi = [];
scale3 = 5/8;

ax1 = subplot(1,2,1);
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax1,colorcet('C06'));
hold on;
vxRaw1a = squeeze(vxRawt(:,:,frame));
vyRaw1a = squeeze(vyRawt(:,:,frame));
vxRaw2a = nan(size(vxRaw1a));
vyRaw2a = nan(size(vyRaw1a));
skip = 6;
zoom_scale = 3;
vxRaw2a(1:skip:end,1:skip:end) = vxRaw1a(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2a(1:skip:end,1:skip:end) = vyRaw1a(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2a(not(BW1)) = nan;
vyRaw2a(not(BW1)) = nan;
imH1Raw4 = quiver(vxRaw2a,vyRaw2a,'k','lineWidth',1,'autoScale','off');
hold on;
axis image; axis off;
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'k';
hold on;
hemi = 'left';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);    
set(gca,'Ydir','reverse')
axis image; axis off;
%
load(fullfile(data_folder,'revision','plane_wave','border_roi2.mat'));
hold on;
% mean flow angle in rectangle roi
vxRaw3 = vxRaw1a(bwroi);
vyRaw3 = vyRaw1a(bwroi);
vxRaw3_mean = mean(vxRaw3);
vyRaw3_mean = mean(vyRaw3);
% mean rectangle center
[row,col] = find(bwroi);
row_mean = mean(row);
col_mean = mean(col);
positions = roi.Position;
positions(end+1,:) = positions(1,:);
plot(positions(:,1),positions(:,2),'r','lineWidth',2);
%%
load(fullfile(data_folder,'revision','plane_wave','angle_mean_all3.mat'));
traceAmp_mean1 = traceAmp_mean;
angle_mean1 = angle(vxy_all);
amp_mean1 = abs(vxy_all);
threshold =0.6;
edges = [-pi:pi/12:pi];
for i = 1:15
    clear angle_mean N
    angle_mean_temp = angle_mean1(i,:);
    angle_mean_filt = angle_mean_temp(amp_mean1(i,:)>threshold);
    [N,edges] =histcounts(angle_mean_filt,edges);
    N1(i,:) = N./size(angle_mean_filt,2);
end
N_mean = mean(N1,1);
N_sem = std(N1,1)./sqrt(15);
subplot(1,2,2);
polarhistogram('BinEdges',edges,'BinCounts',N_mean+N_sem,'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
hold on;
polarhistogram('BinEdges',edges,'BinCounts',N_mean,'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
hold on;
polarhistogram('BinEdges',edges,'BinCounts',N_mean-N_sem,'FaceColor',[0.2,0.2,0.2],'EdgeColor','k');
ax = gca;
ax.ThetaDir = 'clockwise';
%%
print(hs8p, fullfile(save_folder,'FigS8p_border_planar_wave2.pdf'),...
    '-dpdf', '-bestfit', '-painters');
