function getWhiskerMeanMaps(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
BW2 = BW(1:8:end,1:8:end);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
%% right SSp index
clear areaPath
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
hemi = 'right';
Utransformed = projectedAtlas1;
scale = 1;
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp2 = [col,row];
%%
spiral_folder = fullfile(data_folder,'whisker','spirals_all','spirals_grouping');
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
%%
wf_mean_all = [];
for kk = 1:5
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    block_en = T.block_folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    % load data and block
    session_root = fullfile(data_folder,'whisker','task_svd',fname);
    [U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    U = U./mimg;
    bfd(1,1) = T{kk,7};
    bfd(1,2) = T{kk,8};
    expDir = dir(fullfile(session_root,'*_Block.mat'));
    load(fullfile(expDir.folder,expDir.name));
    %
    sigName = 'rewardValve';
    load(fullfile(session_root,[sigName '_raw.mat']));                     % load pd
    load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));     % load tlTimes
    tt = tsToT(tlTimes, numel(pd)); 
    %
    [allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
    flipsUp_iti = diff(flipsUp);
    flipsUp_iti = [0;flipsUp_iti];
    flipsUp = flipsUp(flipsUp<2700);
    %%
    load(fullfile(data_folder,'whisker','rfmap',[fname '.mat']));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
    %%
    for ii = 1:numel(flipsUp)
        stimOn(ii,1) = find(t-flipsUp(ii)>0, 1, 'first');
    end
    frames1 = [-70:70];
    frames  = repmat(frames1,[numel(stimOn),1]);
    frames_stimOn = frames+stimOn;   
    %%
    clear Va indx
    for i = 1:numel(flipsUp)
        ta = t-flipsUp(i);
        indx(i,1) = find(ta>0, 1,'first');
        tt1(i,:) = t(indx(i,1)-70:indx(i,1)+70)-flipsUp(i);
        Va(:,:,i) = dV(:,indx(i,1)-70:indx(i,1)+70);
    end
    Ut= imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
    Ut1 = Ut(1:8:end,1:8:end,:);
    va1 = reshape(Va,size(Va,1),size(Va,2)*size(Va,3));
    wf = reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3))*va1;
    wf = reshape(wf, size(Ut1,1),size(Ut1,2),size(Va,2),numel(flipsUp));  
    if strcmp(T.hemisphere{kk},'left')        
        wf = flip(wf,2);
    end
    %% all trials
    wf_mean =  squeeze(mean(wf,4));
    wf_mean_all(:,:,:,kk) = wf_mean;
    close all;
end
wf_mean2 = mean(wf_mean_all,4);
save(fullfile(save_folder,'whisker_spirals_mean_all.mat'),'wf_mean2','wf_mean_all');
