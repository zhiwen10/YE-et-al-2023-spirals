function h5e = plotSpiralRateAll4(data_folder,save_folder)
%%
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
% freq = [0.5, 2];
freq = [2 8];
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
%%
spiral_count_sum_all =  zeros(3,4,141);
%%
spiral_contra = [];
spiral_ipsi = [];
for i = 1:4
    %%
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder1,'sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    sessions = size(T1,1);
    load(fullfile(data_folder2,[mn '_spirals_task_sort.mat']));
    load(fullfile(data_folder1,'task_outcome',[mn '_task_outcome']));
    for kk = 1:2
    %% left stim
        clear index spiral_cell_temp spiral_cell_temp2
        index = (T_all.label == labels{kk} & T_all.left_contrast > 0);
        spiral_cell_temp = spiral_all(index,:);
        trialN1 = sum(index);
        spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
        for n = 1:141
            clear spiral_temp indx1 indx2
            spiral_temp =  spiral_cell_temp2{n};
            indx1 = (spiral_temp(:,1) > 570 & spiral_temp(:,3)==100);
            spiral_Lstim_contra(n,1) = sum(indx1);
            indx2 = (spiral_temp(:,1) < 570 & spiral_temp(:,3)==100);
            spiral_Lstim_ipsi(n,1) = sum(indx2);
        end
        %% right stim
        clear index spiral_cell_temp spiral_cell_temp2
        index = (T_all.label == labels{kk} & T_all.right_contrast > 0);
        spiral_cell_temp = spiral_all(index ,:);
        trialN2 = sum(index);
        spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
        for n = 1:141
            clear spiral_temp indx1 indx2
            spiral_temp =  spiral_cell_temp2{n};
            indx1 = (spiral_temp(:,1) < 570 & spiral_temp(:,3)==100);
            spiral_Rstim_contra(n,1) = sum(indx1);
            indx2 = (spiral_temp(:,1) > 570 & spiral_temp(:,3)==100);
            spiral_Rstim_ipsi(n,1) = sum(indx2);
        end
        %%
        spiral_contra(:,i,kk) = (spiral_Lstim_contra+spiral_Rstim_contra)./(trialN1+trialN2);
        spiral_ipsi(:,i,kk) = (spiral_Lstim_ipsi+spiral_Rstim_ipsi)./(trialN1+trialN2);
    end
    %%
    kk = 3;
    clear index spiral_cell_temp spiral_cell_temp2
    index = (T_all.label == "miss");
    spiral_cell_temp = spiral_all(index,:);
    trialN = sum(index);
    spiral_cell_temp2 = getConcatTrials(spiral_cell_temp);
    for n = 1:141
        clear spiral_temp indx1 indx2
        spiral_temp =  spiral_cell_temp2{n};
        indx1 = (spiral_temp(:,3)==100);
        spiral_count_sum_miss(n,i,kk) = sum(indx1)/trialN*1;
    end
    
end
%%
spiral_contra = spiral_contra*35;
spiral_ipsi = spiral_ipsi*35;
%%
spiral_count_sum2 = squeeze(sum(spiral_contra(75:82,:,:),1))/8;
spiral_count_sum3 = squeeze(sum(spiral_contra(63:70,:,:),1))/8;
for i = 1:2
    [h1a(i,1),p1a(i,1)] = ttest(spiral_count_sum2(:,i)',spiral_count_sum3(:,i)');
end
%%
spiral_count_sum2 = squeeze(sum(spiral_ipsi(75:82,:,:),1))/8;
spiral_count_sum3 = squeeze(sum(spiral_ipsi(63:70,:,:),1))/8;
for i = 1:2
    [h2a(i,1),p2a(i,1)] = ttest(spiral_count_sum2(:,i)',spiral_count_sum3(:,i)');
end
%%
spiral_count_sum = spiral_ipsi+ spiral_contra;
spiral_count_sum2 = squeeze(sum(spiral_count_sum(75:82,:,:),1))/8;
spiral_count_sum3 = squeeze(sum(spiral_count_sum(63:70,:,:),1))/8;
for i = 1:2
    [h3a(i,1),p3a(i,1)] = ttest(spiral_count_sum2(:,i)',spiral_count_sum3(:,i)');
end
%%
figure;
for i = 1:4
    subplot(4,3,1+3*(i-1));
    plot(1:141,spiral_ipsi(:,i));
    ylim([0,8]);
    hold on;
    %%
    subplot(4,3,2+3*(i-1));
    plot(1:141,spiral_contra(:,i));
    ylim([0,8]);
    hold on;
    %%
    subplot(4,3,3+3*(i-1));
    plot(1:141,spiral_count_sum(:,i));
    ylim([0,8]);
    hold on;
end
%%
spiral_count_sum(:,:,3) = spiral_count_sum_miss(:,:,3)*35;
%%
% spiral_count_sum_all = spiral_count_sum_all*35;
color1 = {'k','g','r'};
h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
for i = 1:3
    spiral_count_sum1 = squeeze(spiral_count_sum(:,:,i));
    spiral_count_mean = squeeze(mean(spiral_count_sum1,2));
    spiral_count_sem = squeeze(std(spiral_count_sum1,[],2))./sqrt(4);
    subplot(1,3,i);
    shadedErrorBar(1:71, spiral_count_mean(35:105), spiral_count_sem(35:105), 'lineprops', color1{i});
    xline(36,'--k');
    xticks([1,18,36,54,71]);
    xticklabels([-1:0.5:1]);
    ylim([0,6]);
    xlim([1,71]);
end
%%
print(h2, 'Task_spirals_pre_post_stim_all_mice_large_spiral2', ...
    '-dpdf', '-bestfit', '-painters');
%%
spiral_count_sum2 = squeeze(sum(spiral_count_sum_all(:,:,75:82),3))/8;
spiral_count_sum3 = squeeze(sum(spiral_count_sum_all(:,:,63:70),3))/8;
[h1,p1] = ttest(spiral_count_sum2',spiral_count_sum3');
color1 = {'k','g','r'};
h2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
for i = 1:3
    subplot(1,3,i);
    scatter(ones(1,4),spiral_count_sum2(i,:),4,'k');
    hold on;
    scatter(ones(1,4)*2,spiral_count_sum3(i,:),4,'k');
    hold on;
    plot([ones(1,4);ones(1,4)*2],[spiral_count_sum2(i,:);spiral_count_sum3(i,:)],'k');
    % ylim([0,3.5]);
    xticks([1,2]);
    xticklabels({'Pre','post'});
    ylabel('Rotating waves/s');
end
%%
print(h2, 'Task_spirals_pre_post_stim_stats', ...
    '-dpdf', '-bestfit', '-painters');