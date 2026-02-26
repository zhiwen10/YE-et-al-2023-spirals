function getBorderPlaneWave3(data_folder,save_folder)
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
%%
epoch = 21000:31500;
useGPU = 1;
load(fullfile(data_folder,'revision','plane_wave','border_roi2.mat'));
%%
vxy_all = [];
traceAmp_mean = [];
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
    traceAmpta = reshape(traceAmpt,size(traceAmpt,1)*size(traceAmpt,2),size(traceAmpt,3));
    traceAmpta = traceAmpta(logical(bwroi(:)),:);
    traceAmp_mean(kk,:) = mean(traceAmpta,1);
    %%
    vxyt = complex(vxRawt,vyRawt);
    vxyt2  = reshape(vxyt,size(vxyt,1)*size(vxyt,2),size(vxyt,3));
    vxyt3 = vxyt2(logical(bwroi(:)),:);  
    % vxyt3 = vxyt3./abs(vxyt3);                                                 % normalize to unit vector
    vxy_all(kk,:) = sum(vxyt3,1)./sum(abs(vxyt3),1);                                        % get sync index
end
%
save(fullfile(save_folder,'angle_mean_all3.mat'),'vxy_all','traceAmp_mean');