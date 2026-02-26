function getCorrectSpiralDensity(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
%%
spiral_folder = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
labels = {"correct","incorrect","miss"};
frames_pre = 62:68; frames_post = 76:82;
framesN = numel(frames_pre);
trialN = zeros(3,1);
spiral_temp_pre_all = cell(3,1);
spiral_temp_post_all = cell(3,1);
for kk = 1:4
    %%
    mn = fnames{kk};
    load(fullfile(spiral_folder,[mn '_spirals_task_sort.mat']),'spiral_all','T_all');
    for j = 1:3
        clear indx spiral_temp spiral_temp1 
        if j ==1
            indx = (T_all.label == labels{j} & (T_all.left_contrast-T_all.right_contrast)>0);
        else
            indx = (T_all.label == labels{j}& abs(T_all.left_contrast-T_all.right_contrast)>0);
        end
        spiral_temp = spiral_all(indx,:);    
        spiral_temp1 = getConcatTrials(spiral_temp);
        spiral_temp_pre = cat(1,spiral_temp1{frames_pre});
        spiral_temp_post = cat(1,spiral_temp1{frames_post});
        %%
        spiral_temp_pre_all{j} = [spiral_temp_pre_all{j};spiral_temp_pre];
        spiral_temp_post_all{j} = [spiral_temp_post_all{j};spiral_temp_post];
        trialN(j,1) = trialN(j,1)+sum(indx);
    end
end

spiral_temp_pre_all2 = cell(3,1);
spiral_temp_post_all2 = cell(3,1);
for j = 1:3
    spiral_temp_pre_all1 = spiral_temp_pre_all{j};
    spiral_temp_pre_all2{j} = spiral_temp_pre_all1(...
        spiral_temp_pre_all1(:,3)==100 &spiral_temp_pre_all1(:,4)>=-2,:);
    %%
    spiral_temp_post_all1 = spiral_temp_post_all{j};
    spiral_temp_post_all2{j} = spiral_temp_post_all1(...
        spiral_temp_post_all1(:,3)==100 &spiral_temp_post_all1(:,4)>=-2,:);
end

hist_bin = 40;
unique_spirals_pre = cell(3,1);
unique_spirals_post = cell(3,1);
for j = 1:3
    unique_spirals_pre{j,1} = density_color_plot2(spiral_temp_pre_all2{j},hist_bin);
    unique_spirals_post{j,1} = density_color_plot2(spiral_temp_post_all2{j},hist_bin);
end
%%
save(fullfile(save_folder,'spiral_density_map_trial_types.mat'),...
    'unique_spirals_pre','unique_spirals_post','trialN');