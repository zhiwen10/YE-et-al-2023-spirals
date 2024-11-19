function [ratio] = getSpiralLargeRatio(T,data_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
for kk = 1:size(T,1)
    clear spiralsT index
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,'spirals\spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    load(fullfile(data_folder,'spirals\rf_tform',...
        [fname '_tform.mat']));
    %%
    t = readNPY(fullfile(data_folder,'spirals\svd',fname,...
        'svdTemporalComponents_corr.timestamps.npy'));                     % read how many total frames
    nframe = numel(t);   
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration >= 2);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));
    filteredSpirals(:,1:2) = round(spiralsT); 
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % only use spirals detected within brain ROI
    filteredSpirals = filteredSpirals(lia,:);
    %% only check spirals >1.7mm
    index = (filteredSpirals(:,3)>=100);
    ratio(kk,1) = sum(index)./nframe;
end
%%
ratio_mean = mean(ratio);
ratio_sem = std(ratio)./sqrt(15);