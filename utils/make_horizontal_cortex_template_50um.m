data_folder = 'E:\spiral_data_share\data';
save_folder = 'E:\spiral_data_share\data\tables';
make_horizontal_cortex_template_50um(data_folder,save_folder);

function make_horizontal_cortex_template_50um(data_folder,save_folder)
%% only select cortex in the atlas
tv = readNPY(fullfile(data_folder,'tables',...
    'template_volume_10um.npy'));                                          % grey-scale "background signal intensity" 
av = readNPY(fullfile(data_folder,'tables',...
    'annotation_volume_10um_by_index.npy'));                               % the number at each pixel labels the area, see note below
%%
template1 =zeros(1320,1140);
for k1 = 1:1320
    for k2 = 1:1140
        secondIndex = find(av(k1,:,k2)>1,60,'first');
        if not(isempty(secondIndex))
            template1(k1,k2) = tv(k1,secondIndex(end),k2);
        end
    end
end
template1 = template1(1:8:end,1:8:end);
save(fullfile(save_folder,'horizontal_cortex_template_50um.mat'),'template1');
end