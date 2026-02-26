function h5e = plotSpiralRateAll3(data_folder,save_folder)
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
% freq = [0.5, 2];
freq = [2 8];
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
%%
spiral_count_sum_all =  zeros(3,4,141);
for i = 1:4
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder1,'sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    sessions = size(T1,1);
    % [spiral_all] = getTaskTrialSpirals(T1); 
    load(fullfile(data_folder2,[mn '_spirals_task_sort.mat']),'spiral_all');
    load(fullfile(data_folder1,'task_outcome',[mn '_task_outcome']));
    %%
    for id = 1:3
        clear index sprial_label
        label = labels{id};   
        index = (T_all.label == label & abs(T_all.left_contrast-T_all.right_contrast)> 0);
        spiral_label = spiral_all(index ,:);
        trialN = sum(index);
        %%
        spiral_count_sum =zeros(141,1);
        for n = 1:141
            clear spiral_temp spiral_temp2 indx1 indx2 indx3 indx 
            spiral_temp = spiral_label(:,n);
            spiral_temp2 = cat(1,spiral_temp{:});
            indx3 = (spiral_temp2(:,3)==100);
            indx = (indx3);
            spiral_count_sum(n,1) = sum(indx)/trialN*1;
        end
        spiral_count_sum_all(id,i,:) = spiral_count_sum;
    end
end
spiral_count_sum_all = spiral_count_sum_all*35;
% spiral_count_sum_all = spiral_count_sum_all./area_cortex;
% spiral_count_sum_all = spiral_count_sum_all./2;
%%
save('task_spiral_count_radius_100','spiral_count_sum_all','labels');
%%
load('task_spiral_count_radius_100');
%%
% spiral_count_sum_all = spiral_count_sum_all*35;
color1 = {'k','g','r'};
h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
for i = 1:3
    spiral_count_sum_all1 = squeeze(spiral_count_sum_all(i,:,:));
    spiral_count_mean = squeeze(mean(spiral_count_sum_all1,1));
    spiral_count_sem = squeeze(std(spiral_count_sum_all1,[],1))./sqrt(4);
    subplot(1,3,i);
    shadedErrorBar(1:71, spiral_count_mean(35:105), spiral_count_sem(35:105), 'lineprops', color1{i});
    xline(36,'--k');
    xticks([1,18,36,54,71]);
    xticklabels([-1:0.5:1]);
    ylim([0,3.5]);
    xlim([1,71]);
end
%%
print(h2, 'Task_spirals_pre_post_stim_all_mice_large_spiral2', ...
    '-dpdf', '-bestfit', '-painters');
%%
spiral_count_sum2 = squeeze(sum(spiral_count_sum_all(:,:,75:82),3))/8;
spiral_count_sum3 = squeeze(sum(spiral_count_sum_all(:,:,63:70),3))/8;
[h1,p1] = ttest(spiral_count_sum2',spiral_count_sum3');
%%
color1 = {'k','g','r'};
h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
for i = 1:3
    subplot(1,3,i);
    scatter(ones(1,4),spiral_count_sum2(i,:),4,'k');
    hold on;
    scatter(ones(1,4)*2,spiral_count_sum3(i,:),4,'k');
    hold on;
    plot([ones(1,4);ones(1,4)*2],[spiral_count_sum2(i,:);spiral_count_sum3(i,:)],'k');
    ylim([0,3.5]);
    xticks([1,2]);
    xticklabels({'Pre','post'});
    ylabel('Rotating waves/s');
end
print(h2, 'Task_spirals_pre_post_stim_stats', ...
    '-dpdf', '-bestfit', '-painters');