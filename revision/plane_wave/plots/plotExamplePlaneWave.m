function hs8n = plotExamplePlaneWave(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
load(fullfile(data_folder,'revision','plane_wave','area_mask.mat'));
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
%
epoch = 990*35:992*35;
frame = 11;
[~,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1(:,:,1:50),dV(1:50,epoch),t,params,freq,rate);
% clear U U1 V dV
tracePhase1 = permute(tracePhase1,[3,1,2]);
useGPU = 1;
[vxRaw,vyRaw] = HS_flowfield(tracePhase1,useGPU);
%% now let's transform flow field to registered atlas
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
hs8n = figure('Renderer', 'painters', 'Position', [50 50 900 700]);
frame = 11;
lineColor = 'k';
hemi = [];
scale3 = 5/8;
BW_all = cat(3,BW_SSp_left,BW_SSp_right,BW_MO_right,BW_SSp_right);
hemi_all = {'left','right','right','right'};
for kk = 1:4
    BW2 = squeeze(BW_all(:,:,kk));
    hemi = hemi_all{kk};
    %%
    ax1 = subplot(2,2,kk);
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
    vxRaw2a(not(BW2)) = nan;
    vyRaw2a(not(BW2)) = nan;
    imH1Raw4 = quiver(vxRaw2a,vyRaw2a,'k','lineWidth',1,'autoScale','off');
    hold on;
    axis image; axis off;
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    scale3 = 5/8;
    lineColor = 'k'; lineColor1 = 'w';
    hold on;
    if ismember(kk,[1,2,4])
        plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
        plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
        plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
        plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
    elseif kk == 3
        plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
    end
    set(gca,'Ydir','reverse')
    axis image; axis off;
end
print(hs8n, fullfile(save_folder,'FigS8n_wave_symmetry_example.pdf'),...
    '-dpdf', '-bestfit', '-painters');