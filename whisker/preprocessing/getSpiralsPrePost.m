function getSpiralsPrePost(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
spiral_folder = fullfile(data_folder,'whisker','spirals_all','spirals_grouping');
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
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
Utransformed = projectedAtlas1;
scale = 1;

hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp_right = [col,row];

hemi = 'left';
[indexSSp2,UselectedSSp2] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexSSp2);
indexSSp_left = [col2,row2];
%%
% areaPath(1) = "/997/8/567/688/695/315/500/985/"; %MOp
areaPath(1) = "/997/8/567/688/695/315/500/993/"; % MOs
MOsArea = strcat(areaPath(:));

% sensoryArea = [];
Utransformed = projectedAtlas1;
scale = 1;
hemi = 'right';
[indexMOs,UselectedSSp] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexMOs);
indexMOs_right = [col,row];

hemi = 'left';
[indexMOs,UselectedSSp2] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexMOs);
indexMOs_left = [col2,row2];
%%
areaPath(1) = "/997/8/567/688/695/315/500/985/"; %MOp
MOsArea = strcat(areaPath(:));

% sensoryArea = [];
Utransformed = projectedAtlas1;
scale = 1;
hemi = 'right';
[indexMOp,UselectedSSp] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexMOp);
indexMOp_right = [col,row];

hemi = 'left';
[indexMOp,UselectedSSp2] = select_area(MOsArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexMOp);
indexMOp_left = [col2,row2];
%%
spirals_pre_all = [];
spirals_post_all = [];
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
    %
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
    %
    load(fullfile(data_folder,'whisker','rfmap',[fname '.mat']));
    %
    for ii = 1:numel(flipsUp)
        stimOn(ii,1) = find(t-flipsUp(ii)>0, 1, 'first');
    end
    frames1 = [-70:70];
    frames  = repmat(frames1,[numel(stimOn),1]);
    frames_stimOn = frames+stimOn;
    % load spirals
    load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));    
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    pwAll = cell2mat(groupedCells); 
    % organize sprials by frame
    [pwAll(:,1),pwAll(:,2)] = transformPointsForward(...
        tform,pwAll(:,1),pwAll(:,2)); 
    pwAll(:,1:2) = round(pwAll(:,1:2));
    if strcmp(T.hemisphere{kk},'left')        
       pwAll(:,1) = 1140-pwAll(:,1);
       pwAll(:,4) = -pwAll(:,4);
    end
    [lia,locb] = ismember(pwAll(:,1:2),brain_index,'rows');
    pwAll = pwAll(lia,:);
    % pwAll = pwAll(pwAll(:,3)>=60,:);
    %%
    spiral_cell = cell(size(frames_stimOn));
    for j = 1:size(frames_stimOn,1)
        for k = 1:size(frames_stimOn,2)
            clear indx1 current_spirals
            current_frame = frames_stimOn(j,k);
            indx1 = find(pwAll(:,5) == current_frame);
            current_spirals = pwAll(indx1,:);
            spiral_cell{j,k} = current_spirals;
        end
    end  
    %%
    spirals_pre = cat(1,spiral_cell{:,67:71});
    spirals_post = cat(1,spiral_cell{:,72:76});
    %%
    spirals_pre_cell{kk,1} = spirals_pre;
    spirals_post_cell{kk,1} = spirals_post;
end
spirals_pre_all = cat(1,spirals_pre_cell{:});
spirals_post_all = cat(1,spirals_post_cell{:});
%%
save(fullfile(save_folder,'Whisker_spirals_pre_post.mat'),...
    'spirals_pre_all','spirals_post_all','spirals_pre_cell','spirals_post_cell');