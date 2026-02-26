function getSpiralsPeriStim(data_folder,save_folder)
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
sensoryArea = [];
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
    %%
    [U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
    %%
    expDir = dir(fullfile(session_root,'*_Block.mat'));
    load(fullfile(expDir.folder,expDir.name));
    sigName = 'rewardValve';
    load(fullfile(session_root,[sigName '_raw.mat']));                     % load pd
    load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));     % load tlTimes
    tt = tsToT(tlTimes, numel(pd)); 
    [allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
    flipsUp_iti = diff(flipsUp);
    flipsUp_iti = [0;flipsUp_iti];
    flipsUp = flipsUp(flipsUp<2700);
    load(fullfile(data_folder,'whisker','rfmap',[fname '.mat']));
    for ii = 1:numel(flipsUp)
        stimOn(ii,1) = find(t-flipsUp(ii)>0, 1, 'first');
    end
    frames1 = [-35:35];
    frames  = repmat(frames1,[numel(stimOn),1]);
    frames_stimOn = frames+stimOn;
    %% load spirals
    load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));    
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    pwAll = cell2mat(groupedCells); 
    %% organize sprials by frame
    [pwAll(:,1),pwAll(:,2)] = transformPointsForward(...
        tform,pwAll(:,1),pwAll(:,2)); 
    pwAll(:,1:2) = round(pwAll(:,1:2));
    [lia,locb] = ismember(pwAll(:,1:2),brain_index,'rows');
    pwAll = pwAll(lia,:);
    %%
    if strcmp(T.hemisphere{kk},'left')        
       pwAll(:,1) = 1140-pwAll(:,1);
       pwAll(:,4) = -pwAll(:,4);
    end
    %%
    spiral_count_sum_left(kk,:,:) = sortSpirals(pwAll,indexSSp_left,frames_stimOn);
    spiral_count_sum_right(kk,:,:) = sortSpirals(pwAll,indexSSp_right,frames_stimOn);
end
spiral_count_sum_left = spiral_count_sum_left*35;
spiral_count_sum_right = spiral_count_sum_right*35;
%%
save(fullfile(save_folder,'Whisker_spirals_peri_stimulus.mat'),...
    'spiral_count_sum_left','spiral_count_sum_right');
