function h5e = plotSpiralRateAll(data_folder,save_folder)
%%
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
%%
sampleN = 71;
half_spiral_rate_all = zeros(4,3,sampleN);
for i = 1:4
    %%
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    %%
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome.mat']));
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_trial_ID.mat']));
    load(fullfile(data_folder,'task','spirals_half',[mn '_half_spirals4.mat']));
    for j = 1:numel(labels)
        if j==1
            indx = (T_all.label== labels{j}& (T_all.left_contrast-T_all.right_contrast)>=0.25);
        else
            indx = (T_all.label== labels{j});
        end
        halfSpiral_temp = halfSpiral_all(indx,:,:);
        %%
        indx2 = (halfSpiral_temp(:,:,2)==-1);
        cw_sum = sum(indx2,1);
        cw_rate = cw_sum./sum(indx);
        half_spiral_rate_all(i,j,:) = cw_rate;
    end
end
half_rate_mean = squeeze(mean(half_spiral_rate_all, 1));
half_rate_sem = squeeze(std(half_spiral_rate_all,[],1))./sqrt(4);
% half spiral stats
pre_rate = squeeze(half_spiral_rate_all(:,:,36));
post_rate = squeeze(half_spiral_rate_all(:,:,45));
[h2,p2] = ttest(pre_rate,post_rate);
%% full spiral raster plot
load(fullfile(data_folder,'task','spirals_large','task_spiral_count_radius_100.mat'));
% full spiral stats
spiral_count_mean2 = squeeze(mean(spiral_count_sum_all(:,:,75:82),3));
spiral_count_mean3 = squeeze(mean(spiral_count_sum_all(:,:,63:70),3));
[h3,p3] = ttest(spiral_count_mean2',spiral_count_mean3');

spiral_count_mean_pre1 = squeeze(mean(spiral_count_sum_all(1,:,63:70),3));
spiral_count_mean_post1 = squeeze(mean(spiral_count_sum_all(1,:,75:82),3));
spira_mean_all_pre1 = mean(spiral_count_mean_pre1,2);
spira_mean_all_post1 = mean(spiral_count_mean_post1,2);
spiral_sem_pre1 = squeeze(std(spiral_count_mean_pre1))./sqrt(4);
spiral_sem_post1 = squeeze(std(spiral_count_mean_post1))./sqrt(4);
%% half spiral raster plot
color1 = {'k','g','r'};
h5e = figure('Renderer', 'painters', 'Position', [100 100 800 500]);
for i = 1:3
    subplot(2,3,i);
    shadedErrorBar(1:71,half_rate_mean(i,:),half_rate_sem(i,:),'lineprops', color1{i});
    ylim([0,0.15]);
    xline(36,'--k');
    xticks([1,18,36,54,71]);
    xticklabels({'-1','-0.5','0','0.5','1'});
    xlim([1,71]);
    yticks([0,0.05,0.1,0.15]);
    yticklabels({'0','5%','10%','15%'});
end

for i = 1:3
    spiral_count_sum_all1 = squeeze(spiral_count_sum_all(i,:,:));
    spiral_count_mean = squeeze(mean(spiral_count_sum_all1,1));
    spiral_count_sem = squeeze(std(spiral_count_sum_all1,[],1))./sqrt(4);
    subplot(2,3,i+3);
    shadedErrorBar(1:71, spiral_count_mean(35:105), spiral_count_sem(35:105), 'lineprops', color1{i});
    xline(36,'--k');
    xticks([1,18,36,54,71]);
    xticklabels({'-1','-0.5','0','0.5','1'});
    ylim([0,0.06]);
    xlim([1,71]);
    yticks([0,0.02,0.04,0.06]);
    yticklabels({'0','2%','4%','6%'});
end

%%
print(h5e, fullfile(save_folder,'Fig5e_half_and_full_spirals'),...
    '-dpdf', '-bestfit', '-painters');
