function getBorderPlaneWave(data_folder,save_folder)
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
% use 600s to 635s as epoch examples, 600x35 = 21000 samples
epoch = 21000:31500;
% epoch = 21000:22050;
useGPU = 1;
% load(fullfile(data_folder,'revision','planar_wave','SSP_MOp_boarder_roi.mat'));
load(fullfile(data_folder,'revision','plane_wave','RSP_roi.mat'));
bwroi_RSP = bwroi;
load(fullfile(data_folder,'revision','plane_wave','ACC_roi.mat'));
bwroi_ACC = bwroi;
load(fullfile(data_folder,'revision','plane_wave','VISp_roi.mat'));
bwroi_VISp = bwroi;
load(fullfile(data_folder,'revision','plane_wave','SSP_MOp_border_roi.mat'));
bwroi_border = bwroi;
%%
for kk = 1:15
    tic
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
    %% flow field 
    % get flow field from unregistered frame frist, then transform to
    % registered, to avoid interpolation problem with circular phase 
    U1 = U(1:params.downscale:end,1:params.downscale:end,:);
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    U1 = U1./mimg1;
    [~,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1(:,:,1:50),dV(1:50,epoch),t,params,freq,rate);
    tracePhase1 = permute(tracePhase1,[3,1,2]);
    % tracePhase2 = tracePhase1(epoch,:,:);
    [vxRaw,vyRaw] = HS_flowfield(tracePhase1,useGPU);
    % now let's transform flow field to registered atlas
    BW1 = BW(1:params.downscale:end,1:params.downscale:end);
    vxRaw1 = permute(vxRaw,[2,3,1]);
    vyRaw1 = permute(vyRaw,[2,3,1]);
    vxRawt = imwarp(vxRaw1,tform,'OutputView',imref2d(size(BW1)));
    vyRawt = imwarp(vyRaw1,tform,'OutputView',imref2d(size(BW1)));
    traceAmpt = imwarp(traceAmp1,tform,'OutputView',imref2d(size(BW1)));
    %%
    [angle_mean_border(kk,:),traceAmp_mean_border(kk,:)] = getFlowAmpROI(traceAmpt,vxRawt,vyRawt,bwroi_border);
    [angle_mean_RSP(kk,:),traceAmp_mean_RSP(kk,:)] = getFlowAmpROI(traceAmpt,vxRawt,vyRawt,bwroi_RSP);
    [angle_mean_ACC(kk,:),traceAmp_mean_ACC(kk,:)] = getFlowAmpROI(traceAmpt,vxRawt,vyRawt,bwroi_ACC);
    [angle_mean_VISp(kk,:),traceAmp_mean_VISp(kk,:)] = getFlowAmpROI(traceAmpt,vxRawt,vyRawt,bwroi_VISp);
end
%%
save(fullfile(save_folder,'angle_mean_all.mat'),'angle_mean_border','traceAmp_mean_border',...
    'angle_mean_RSP','traceAmp_mean_RSP','angle_mean_ACC','traceAmp_mean_ACC',...
    'angle_mean_VISp','traceAmp_mean_VISp');