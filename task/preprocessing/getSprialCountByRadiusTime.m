function getSprialCountByRadiusTime(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
pix_sum = sum(BW(:));
area_cortex = pix_sum*pixArea;
%%
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
radius = 40:10:100;
spiral_count_sum_all =  zeros(3,4,7,141);
%%
for id = 1:3
    %%
    label = labels{id};
    %%
    for i = 1:4
        %%
        mn = fnames{i};
        T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
        T1 = T_session(T_session.label == "task",:);
        T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
        sessions = size(T1,1);
        load(fullfile(data_folder,'task','spirals',[mn '_spirals_task_sort.mat']));
        load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
        %%    
        index = (T_all.label == label & abs(T_all.left_contrast-T_all.right_contrast)>0);
        spiral_label = spiral_all(index,:);
        trialN = sum(index);
        %%
        for kk = 1:7
            spiral_count =zeros(size(spiral_label,1),141);
            for m = 1:size(spiral_label,1)
                spiral_correct_temp = spiral_label(m,:);
                for n = 1:141
                    spiral_temp = spiral_correct_temp{n};
                    if not(isempty(spiral_temp))
                        indx = (spiral_temp(:,3)== radius(kk));
                        spiral_temp2 = spiral_temp(indx,:);
                        if not(isempty(spiral_temp2))
                            spiral_count(m,n) = 1;
                        end            
                    end
                end
            end
            spiral_count_sum = sum(spiral_count,1)./trialN*35;
            spiral_count_sum_all(id,i,kk,:) = spiral_count_sum;
        end    
    end
end
%%
save(fullfile(save_folder,'task_spiral_count_by_radius_time'),...
    'spiral_count_sum_all','labels');