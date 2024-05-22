data_folder = 'E:\spiral_data_share\data';
save_folder = 'E:\spiral_data_share\data\tables';
make_horizontal_cortex_atlas_50um(data_folder,save_folder);
%%
function make_horizontal_cortex_atlas_50um(data_folder,save_folder)
[atlas, metaAVGT] = nrrdread(fullfile(data_folder,...
    'tables', 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
save(fullfile(save_folder,'horizontal_cortex_atlas_50um.mat'),'atlas1');
end