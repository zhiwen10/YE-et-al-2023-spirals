function hr1a = plotExamplePlaneWaveSeries(data_folder,save_folder)
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
[U,V,t,mimg] = loadUVt1(session_root); 
dV = [zeros(size(V,1),1) diff(V,[],2)];
% registration
load(fullfile(data_folder,'spirals\rf_tform_8x',[fname '_tform_8x.mat']));   % load atlas transformation matrix tform;
% get flow field from unregistered frame frist, then transform to
% registered, to avoid interpolation problem with circular phase 
U1 = U(1:params.downscale:end,1:params.downscale:end,:);
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
U1 = U1./mimg1;
%%
% epoch = 820*35:822*35;
% kk = 5;
epoch = 990*35:992*35;
[~,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1(:,:,1:50),dV(1:50,epoch),t,params,freq,rate);
trace = squeeze(U1(57,49,1:50))'*dV(1:50,epoch); % right V1
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
vxyt = complex(vxRawt,vyRawt);
vxyt2  = reshape(vxyt,size(vxyt,1)*size(vxyt,2),size(vxyt,3));
vxyta = vxyt2(logical(BW1(:)),:);
vxyta = vxyta(not(abs(vxyta(:,1))==0),:);
vxyta = vxyta./abs(vxyta);
vxytb = sum(vxyta,1)./size(vxyta,1);
vxy_sync = abs(vxytb);
%%
traceAmpt2t = reshape(traceAmp1t,size(traceAmp1t,1)*size(traceAmp1t,2),size(traceAmp1t,3));
traceAmpt3t = mean(traceAmpt2t(logical(BW1(:)),:),1,'omitnan');
traceAmpt3t = traceAmpt3t(1:end-1);
%%
trace = trace(1:end-1);
figure;
plot(0:1/35:(numel(traceAmpt3t)-1)/35,traceAmpt3t*100+1,'r');
hold on;
plot(0:1/35:(numel(traceAmpt3t)-1)/35,trace*100+0.5,'b');
hold on;
plot(0:1/35:(numel(traceAmpt3t)-1)/35,vxy_sync,'k');
%%
hr1a = figure('Renderer', 'painters', 'Position', [50 50 900 400]);
frame = 5;
lineColor = 'k';
hemi = [];
scale3 = 5/8;
a = 40; b = 90; % height
c = 20; d = 70; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
for i = 1:10
    ax1 = subplottight(3,10,i);
    im1 = imagesc(squeeze(tracePhase1t(frame+i,:,:)));
    colormap(ax1,colorcet('C06'));
    hold on;
    axis image; axis off;
    set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
%     hold on;
%     overlayOutlines(coords,params.downscale);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    %%
    ax1 = subplottight(3,10,10+i);
    im1 = imagesc(squeeze(tracePhase1t(frame+i,:,:)));
    colormap(ax1,colorcet('C06'));
    hold on;
    vxRaw1a = squeeze(vxRawt(:,:,frame+i));
    vyRaw1a = squeeze(vyRawt(:,:,frame+i));
    vxRaw2a = nan(size(vxRaw1a));
    vyRaw2a = nan(size(vyRaw1a));
    skip = 5;
    zoom_scale = 3;
    vxRaw2a(1:skip:end,1:skip:end) = vxRaw1a(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2a(1:skip:end,1:skip:end) = vyRaw1a(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2a(not(BW1)) = nan;
    vyRaw2a(not(BW1)) = nan;
    imH1Raw4 = quiver(vxRaw2a,vyRaw2a,'k','lineWidth',1,'autoScale','off');
    hold on;
    axis image; axis off;
    set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
%     hold on;
%     overlayOutlines(coords,params.downscale);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    %%
    ax4(i) = subplottight(3,10,20+i);
    im1 = imagesc(squeeze(tracePhase1t(frame+i,:,:)));
    colormap(ax4(i),colorcet('C06'));
    hold on;
    vxRaw1a = squeeze(vxRawt(:,:,frame+i));
    vyRaw1a = squeeze(vyRawt(:,:,frame+i));
    vxRaw2a = nan(size(vxRaw1a));
    vyRaw2a = nan(size(vyRaw1a));
    skip = 5;
    zoom_scale = 3;
    vxRaw2a(1:skip:end,1:skip:end) = vxRaw1a(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2a(1:skip:end,1:skip:end) = vyRaw1a(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2a(not(BW1)) = nan;
    vyRaw2a(not(BW1)) = nan;
    imH1Raw4 = quiver(vxRaw2a,vyRaw2a,'k','lineWidth',1,'autoScale','off');
    hold on;
    axis image; axis off;
    set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
%     hold on;
%     overlayOutlines(coords,params.downscale);
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    ylim([a,b]);
    xlim([c,d]);
    ax4(i).Position(1) = ax4(i).Position(1)+0.005;
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.01;
    ax4(i).Position(4) = ax4(i).Position(4)-0.025;
end
%%
print(hr1a, fullfile(save_folder,'FigR1a_example_flow_field.pdf'),...
    '-dpdf', '-bestfit', '-painters');