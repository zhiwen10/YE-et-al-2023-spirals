function getPassiveSpiralPrePost(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
frames_pre = 59:68; frames_post = 75:84;
framesN = numel(frames_pre);
high_trialN_all = 0;
spirals_high_pre_all = [];
spirals_high_post_all = [];
for kk = 1:4
    %%
    mn = fnames{kk};
    load(fullfile(data_folder,'task','spirals',[mn '_spirals_passive_sort.mat']));
    high_trialN = spiral_high_stimL_trialN + spiral_high_stimR_trialN;
    %%
    [spiral_high_L] = getConcatTrials(spiral_high_stimL);
    [spiral_high_R] = getConcatTrials(spiral_high_stimR);
    spiral_high_all = combineSpiralsLR(spiral_high_L, spiral_high_R);
    %%
    spirals_high_pre = cat(1,spiral_high_all{frames_pre});
    spirals_high_post = cat(1,spiral_high_all{frames_post});
    spirals_high_pre_all = [spirals_high_pre_all;spirals_high_pre];
    spirals_high_post_all = [spirals_high_post_all;spirals_high_post];
    high_trialN_all = high_trialN_all+high_trialN;
end
%%
hist_bin = 40;
[unique_spirals_high_pre] = density_color_plot2(spirals_high_pre_all,hist_bin);
[unique_spirals_high_post] = density_color_plot2(spirals_high_post_all,hist_bin);
unique_spirals_all = {unique_spirals_high_pre,unique_spirals_high_post};
labels = {'High pre-stim','High post-stim'};
trialN = [high_trialN_all,high_trialN_all];
%%
save(fullfile(save_folder,'PassiveSpiralsPrePost.mat'),'unique_spirals_all','trialN');