function getSpiralSpeedConcat(T, data_folder, save_folder)
grid_x = [250,350];
grid_y = [350,550];
%%
angular_all = [];
distance_all = [];
count1 = 1;
for radiusi = 40:10:100
    for kk = 1:15
        clear filteredSpirals2 angle_offset_all distance_offset_all  
        %%
        mn = T.MouseID{kk};
        tda = T.date(kk);
        en = T.folder(kk);    
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        fname = [mn '_' tdb '_' num2str(en)];
        load(fullfile(data_folder,'spirals_speed',[fname '.mat']));
        load(fullfile(data_folder,'spirals_grouped',...
            [fname '_spirals_group_fftn.mat']));
        %%
        clear index1 indx2 groupedCells
        indx2 = cellfun(@(x) size(x,1), archiveCell);
        groupedCells = archiveCell(indx2>=2);
        filteredSpirals2 = cell2mat(groupedCells);
        %%
        index1 = (filteredSpirals2(:,2)>grid_x(1)...
            & filteredSpirals2(:,2)<grid_x(2)...
            &filteredSpirals2(:,1)>grid_y(1)...
            & filteredSpirals2(:,1)<grid_y(2));
        filteredSpirals2 = filteredSpirals2(index1,:);
        indx2 = (filteredSpirals2(:,3) == radiusi);
        %%
        angle_offset_all1 = cellfun(@(x) mean(x,'omitnan'),...
            angle_offset_all,'UniformOutput',false);
        angle_offset_all_select = angle_offset_all1(indx2);
        angle_offset{kk,count1} = angle_offset_all_select;    
        %%
        distance_offset_all_select = distance_offset_all(indx2);           % radius of 8*10=80 pixels
        distance_offset{kk,count1} = distance_offset_all_select;
    end
    %% 
    count1 = count1+1;
end
%%
save(fullfile(save_folder,'speed_all.mat'),'angle_offset','distance_offset');