function getINCorrectSpiralPrePost2(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
% frames_pre = 59:68; frames_post = 75:84;
% frames_pre = 64:66; frames_post = 78:80;
frames_pre = 63:67; frames_post = 77:81;
incorrect_trialN_all = 0;
spirals_incorrect_pre_all = [];
spirals_incorrect_post_all = [];
%%
for kk = 1:4
    %%
    mn = fnames{kk};
    load(fullfile(data_folder,'task','spirals',[mn '_spirals_task_sort.mat']));
    %%
    incorrect_trialN = incorrect_L_trialN + incorrect_R_trialN;  
    spiral_incorrect_all = combineSpiralsLR(spiral_incorrect_L, spiral_incorrect_R);
    %
    framesN = numel(frames_pre);
    spirals_incorrect_pre = cat(1,spiral_incorrect_all{frames_pre});
    spirals_incorrect_post = cat(1,spiral_incorrect_all{frames_post});
    %
    spirals_incorrect_pre_all = [spirals_incorrect_pre_all;spirals_incorrect_pre];
    spirals_incorrect_post_all = [spirals_incorrect_post_all;spirals_incorrect_post];
    incorrect_trialN_all = incorrect_trialN_all+incorrect_trialN;
end
%%
hist_bin = 40;
[unique_spirals_incorrect_pre] = density_color_plot2(spirals_incorrect_pre_all,hist_bin);
[unique_spirals_incorrect_post] = density_color_plot2(spirals_incorrect_post_all,hist_bin);
unique_spirals_all = {unique_spirals_incorrect_pre,unique_spirals_incorrect_post};
labels = {'Incorrect pre-stim','Incorrect post-stim'};
trialN = [incorrect_trialN_all,incorrect_trialN_all];
%%
save(fullfile(save_folder,'IncorrectSpiralsPrePost2.mat'),'unique_spirals_all','trialN');