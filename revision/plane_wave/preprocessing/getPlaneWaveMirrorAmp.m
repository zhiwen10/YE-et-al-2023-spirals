function getPlaneWaveMirrorAmp(data_folder,save_folder)
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
%%
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
% use 600s to 635s as epoch examples, 600x35 = 21000 samples
epoch = 21000:31500;
% epoch = 21000:22050;
useGPU = 1;
%%
for kk = 1:15
    %%
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
    % now let's transform flow field to registered atlas
    BW1 = BW(1:params.downscale:end,1:params.downscale:end);
    traceAmpt = imwarp(traceAmp1,tform,'OutputView',imref2d(size(BW1)));
    traceAmpt2 = reshape(traceAmpt,size(traceAmpt,1)*size(traceAmpt,2),size(traceAmpt,3));
    %%
    traceAmpt_MO_left(kk,:) = mean(traceAmpt2(logical(BW_MO_left(:)),:),1,'omitnan');
    traceAmpt_MO_right(kk,:) = mean(traceAmpt2(logical(BW_MO_right(:)),:),1,'omitnan');
    traceAmpt_SSp_left(kk,:) = mean(traceAmpt2(logical(BW_SSp_left(:)),:),1,'omitnan');
    traceAmpt_SSp_right(kk,:) = mean(traceAmpt2(logical(BW_SSp_right(:)),:),1,'omitnan');
end
%%
save(fullfile(save_folder,'flow_mirror_amp_all.mat'),'traceAmpt_MO_left','traceAmpt_MO_right',...
    'traceAmpt_SSp_left','traceAmpt_SSp_right');