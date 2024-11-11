function getPlaneWaveMirror_newsession(data_folder,save_folder)
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
load(fullfile(data_folder,'revision','plane_wave','area_mask.mat'));
% use 600s to 635s as epoch examples, 600x35 = 21000 samples
epoch = 21000:31500;
% epoch = 21000:22050;
useGPU = 1;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
for kk = 1:4
    %%
    clear spiralsT filteredSpirals unique_spirals unique_spirals_unit indx2
    mn = fnames{kk};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T = T_session(T_session.label == "passive",:);
    %%
    if kk <= 2
        session = 4;
    else
        session = 5;
    end
    mn = T.MouseID{session};
    tda = T.date(session);
    en = T.folder(session);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'task','task_svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root); 
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    % registration
    load(fullfile(data_folder,'revision','spirals_new','rfmap_8x',[fname '.mat']));   % load atlas transformation matrix tform;
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
    vxyt = complex(vxRawt,vyRawt);
    vxyt2  = reshape(vxyt,size(vxyt,1)*size(vxyt,2),size(vxyt,3));
    %%
    vxyt_MO_left = vxyt2(logical(BW_MO_left(:)),:);
    vxyt_MO_right = vxyt2(logical(BW_MO_right(:)),:);
    vxyt_SSp_left = vxyt2(logical(BW_SSp_left(:)),:);
    vxyt_SSp_right = vxyt2(logical(BW_SSp_right(:)),:);
    %%
    vxyt_MO_left = vxyt_MO_left(not(abs(vxyt_MO_left(:,1))==0),:);
    vxyt_MO_right = vxyt_MO_right(not(abs(vxyt_MO_right(:,1))==0),:);
    vxyt_SSp_left = vxyt_SSp_left(not(abs(vxyt_SSp_left(:,1))==0),:);
    vxyt_SSp_right = vxyt_SSp_right(not(abs(vxyt_SSp_right(:,1))==0),:);
    %%
    vxyt_MO_left_abs = abs(vxyt_MO_left);
    vxyt_MO_right_abs = abs(vxyt_MO_right);
    vxyt_SSp_left_abs = abs(vxyt_SSp_left);
    vxyt_SSp_right_abs = abs(vxyt_SSp_right);
    %%
    vxy_MO_left(kk,:) = sum(vxyt_MO_left,1)./sum(vxyt_MO_left_abs,1);
    vxy_MO_right(kk,:) = sum(vxyt_MO_right,1)./sum(vxyt_MO_right_abs,1);
    vxy_SSp_left(kk,:) = sum(vxyt_SSp_left,1)./sum(vxyt_SSp_left_abs,1);
    vxy_SSp_right(kk,:) = sum(vxyt_SSp_right,1)./sum(vxyt_SSp_right_abs,1);
end
%
save(fullfile(save_folder,'flow_mirror_all_newsession2.mat'),'vxy_MO_left','vxy_MO_right',...
    'vxy_SSp_left','vxy_SSp_right');