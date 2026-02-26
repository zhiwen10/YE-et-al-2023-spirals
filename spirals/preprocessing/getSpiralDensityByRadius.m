function getSpiralDensityByRadius(T,data_folder,save_folder)
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
spirals_all = [];
for kk = 1:size(T,1)
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    load(fullfile(data_folder,'spirals','rf_tform',...
        [fname '_tform.mat']));
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=1);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));
    filteredSpirals(:,1:2) = round(spiralsT); 
    spirals_all = [spirals_all;filteredSpirals];
end
spirals_all(:,1:2) = round(spirals_all(:,1:2));
[lia,locb] = ismember(spirals_all(:,1:2),brain_index,'rows');
spirals_all = spirals_all(lia,:);
%%
radius = 40:10:100;
hist_bin = 40;
for i  =1:numel(radius)
    clear spirals_temp
    radius_i = radius(i);
    spirals_temp = spirals_all(spirals_all(:,3)==radius_i,:);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
    filename = ['histogram_' num2str(radius_i) 'radius.mat'];
    save(fullfile(save_folder,filename), 'unique_spirals');
end