function getSpeedForRadius(T,radius,data_folder,save_folder)
scale1 = 8;
grid_x = [250,350];
grid_y = [350,550];
%%
angular_velocity_all = [];
linear_velocity_all = [];
%%
color_stairs = parula(14);
count1 = 1;
for kk = 1:15
    clear filteredSpirals2 angle_offset_all distance_offset_all ...
        index1 indx2 groupedCells
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,'spirals','spirals_speed',[fname '.mat']));
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    indx2 = cellfun(@(x) size(x,1), archiveCell);
    groupedCells = archiveCell(indx2>=2);
    filteredSpirals2 = cell2mat(groupedCells);
    %%
    index1 = (filteredSpirals2(:,2)>grid_x(1) ...
        & filteredSpirals2(:,2)<grid_x(2)...
        &filteredSpirals2(:,1)>grid_y(1)...
        & filteredSpirals2(:,1)<grid_y(2));
    filteredSpirals2 = filteredSpirals2(index1,:);
    indx2 = (filteredSpirals2(:,3) == radius);
    %%
    angle_offset_all1 = cellfun(@(x) mean(x,'omitnan'),...
        angle_offset_all,'UniformOutput',false);
    angle_offset_all_select = angle_offset_all1(indx2);
    angle_offset_all_sessions{kk} = angle_offset_all_select;    
    distance_offset_all_select = distance_offset_all(indx2);               % radius of 8*10=80 pixels
    distance_offset_all_sessions{kk} = distance_offset_all_select;
end
%%
angle_offset_all_sessions2 = horzcat(angle_offset_all_sessions{:});
distance_offset_all_sessions2 = horzcat(distance_offset_all_sessions{:});
angular_velocity_all = vertcat(angle_offset_all_sessions2{:})*35;
distance_offset_all1 = cellfun(@(x) mean(x,1,'omitnan'),...
    distance_offset_all_sessions2,'UniformOutput',false);
linear_velocity_all = vertcat(distance_offset_all1{:});
%%
save(fullfile(save_folder,['speed_' num2str(radius) 'pixels.mat']),...
    'angular_velocity_all','linear_velocity_all');