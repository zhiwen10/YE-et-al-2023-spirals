%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%%
hs = figure('Renderer', 'painters', 'Position', [100 100 300 250]);
labels = {'hemi_predict','MO_predict'};
for j = 1:2
    %%
    clearvars -except labels j data_folder T hs
    label = labels{j};
    load(fullfile(data_folder,'revision2',label,'spirals_compare',...
        'spiral_compare_sessions_neighbor.mat'));
    match_left = [];
    match_right = [];
    for kk = 1:size(T,1)
        clear match1 match2 match_temp
        left_temp = spiral_left_match_all{kk};
        match1 = left_temp(:,11);
        match_left(kk,1) = numel(match1);
        match_left(kk,2) = sum(match1);
        match_left(kk,3) = sum(match1)/numel(match1);
        %%
        right_temp = spiral_right_match_all{kk};
        match2 = right_temp(:,11);
        match_right(kk,1) = numel(match2);
        match_right(kk,2) = sum(match2);
        match_right(kk,3) = sum(match2)/numel(match2);
    end
    %%
    load(fullfile(data_folder,'revision2',label,'spirals_compare',...
        'spiral_compare_sessions_neighbor_perm.mat'));
    match_left_permute = [];
    match_right_permute = [];
    for kk = 1:size(T,1)
        clear match1 match2 match_temp
        left_temp = spiral_left_match_all_perm{kk};
        match1 = left_temp(:,11);
        match_left_permute(kk,1) = numel(match1);
        match_left_permute(kk,2) = sum(match1);
        match_left_permute(kk,3) = sum(match1)/numel(match1);
        %%
        right_temp = spiral_right_match_all_perm{kk};
        match2 = right_temp(:,11);
        match_right_permute(kk,1) = numel(match2);
        match_right_permute(kk,2) = sum(match2);
        match_right_permute(kk,3) = sum(match2)/numel(match2);
    end
    %%
    subplot(1,2,j);
    scatter(ones(size(match_left(:,3))),match_left(:,3),18,'r','filled');
    hold on;
    scatter(2*ones(size(match_left_permute(:,3))),match_left_permute(:,3),18,'k','filled');
    hold on;
    plot([ones(size(match_left(:,3))),2*ones(size(match_left_permute(:,3)))]',[match_left(:,3),match_left_permute(:,3)]','k');
    hold on;
    yline(0.5,'--k');
    ylim([0.2,1]);
    [ha2,pa2] = ttest(match_left(:,3),match_left_permute(:,3));
    mean_all = mean(match_left(:,3));
    sem_all = std(match_left(:,3))./sqrt(size(match_left_permute,1));
    mean_perm = mean(match_left_permute(:,3));
    sem_perm = std(match_left_permute(:,3))./sqrt(size(match_left_permute,1));
    xlabel('Spiral counts');
    ylabel('Spiral matching ratio');
    text(1.5,0.7,['p =' num2str(pa2)]);
    yticks([0.2:0.2:1])
    yticklabels({'0.2','0.4','0.6','0.8','1.0'})
    xlim([0.8,2.2]);
end
%%
print(hs,'spirals_matching_ratio_all2.pdf','-dpdf', '-bestfit', '-painters');