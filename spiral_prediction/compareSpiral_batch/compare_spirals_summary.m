% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
% codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow_new1';
% T1 = readtable(fullfile(codefolder, 'final_table.csv'));

codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T1 = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
raw_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_raw_fftn';
reg_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\registration';
amp_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\amplitude';
%%
count1 = 1;
session_n = [];
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    clear indx T2
    indx = contains(T1.Area,area(1,kk));
    session_n(count1) = sum(indx);
    count1= count1+1;
end
%%
load('spiral_compare_sessions_neighbor.mat');
load('spiral_compare_sessions_neighbor_permute.mat');
% spiral_left_match_all(22:23) = [];
% spiral_right_match_all(22:23) = [];
% spiral_left_match_all_perm(22:23) = [];
% spiral_right_match_all_perm(22:23) = [];
match_all = [];
match_perm = [];
for kk = 1:29
    clear match1 match2 match_temp
    left_temp = spiral_left_match_all{kk};
    match1 = left_temp(:,11);
    right_temp = spiral_right_match_all{kk};
    match2 = right_temp(:,11);
    match_temp = [match1;match2];   
    match_all(kk,1) = numel(match_temp);
    match_all(kk,2) = sum(match_temp)/numel(match_temp);
    
    clear match1p match2p match_tempp
    left_temp = spiral_left_match_all_perm{kk};
    match1p = left_temp(:,11);
    right_temp = spiral_right_match_all_perm{kk};
    match2p = right_temp(:,11);
    match_tempp = [match1p;match2p];   
    match_perm(kk,1) = numel(match_tempp);
    match_perm(kk,2) = sum(match_tempp)/numel(match_tempp);
end
%%
color2 = {'g','r','c','m'};
ratio_all = {};
ratio_perm = {};
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
count1 = 1;
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    clear indx T2
    indx = contains(T1.Area,area(1,kk));
    T2 = T1(indx,:);
    match_all_temp = match_all(indx,:);
    match_perm_temp = match_perm(indx,:);
    ratio_all{count1} = match_all_temp;
    ratio_perm{count1} = match_perm_temp;
    subplot(1,3,count1);
    scatter(1,match_all_temp(:,2),18,color2{kk},'filled');
    hold on;
    scatter(2,match_perm_temp(:,2),18,'k','filled');
    hold on;
    plot([ones(size(match_all_temp(:,2))),2*ones(size(match_perm_temp(:,2)))]',[match_all_temp(:,2),match_perm_temp(:,2)]','k');
    hold on;
    yline(0.5,'--k');
    ylim([0.3,0.9]);
    [ha2(count1),pa2(count1)] = ttest(match_all_temp(:,2),match_perm_temp(:,2));
    xlabel('Spiral counts');
    ylabel('Spiral matching ratio');
    text(1.5,0.7,['p =' num2str(pa2(count1))]);
    yticks([0.3:0.1:0.9])
    yticklabels({'0.3','0.4','0.5','0.6','0.7','0.8','0.9'})
    
    count1 = count1+1;
end
%%
print(h1, 'sprial_matching_ratio_area', '-dpdf', '-bestfit', '-painters');
%%
mean_area = cellfun(@mean,ratio_all,'UniformOutput',false);
mean_area = cat(1,mean_area{:});
std_area = cellfun(@std,ratio_all,'UniformOutput',false);
std_area = cat(1,std_area{:});
count_area = cellfun(@(x) size(x,1),ratio_all);
count_area = repmat(count_area',1,2);
sem_area = std_area./sqrt(count_area);

mean_area_perm = cellfun(@mean,ratio_perm,'UniformOutput',false);
mean_area_perm = cat(1,mean_area_perm{:});
std_area_perm = cellfun(@std,ratio_perm,'UniformOutput',false);
std_area_perm = cat(1,std_area_perm{:});
count_area = cellfun(@(x) size(x,1),ratio_perm);
count_area = repmat(count_area',1,2);
sem_area_perm = std_area_perm./sqrt(count_area);
%%
for i =1:3
    clear ratio_i ratio_i1
    ratio_i = ratio_all{i};
    ratio_i1 = (ratio_i(:,2)-0.5);
    [hh(i),pp(i)] = ttest(ratio_i1);
end
%% version2
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
scatter(match_all(:,1),match_all(:,2),'r','filled');
hold on;
scatter(match_perm(:,1),match_perm(:,2),'k','filled');
hold on;
plot([match_all(:,1),match_perm(:,1)]',[match_all(:,2),match_perm(:,2)]','k');
ylim([0.3,0.8]);
[ha2,pa2] = ttest(match_all(:,2),match_perm(:,2));
[ha,pa] = ttest(match_perm(:,2)-0.5);
xlabel('Spiral counts');
ylabel('Spiral matching ratio');
%%
color2 = {'g','r','c','m'};
count1 = 1;
% edges = [0:0.00125:0.025];
edges = [0:0.0025:0.025];
edges1 = edges(1,2:end);
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    subplot(1,3,count1);
    clear indx T2
    indx = contains(T1.Area,area(1,kk));
    T2 = T1(indx,:);
    spiral_left_match_all1 = spiral_left_match_all(indx);
    spiral_right_match_all1 = spiral_right_match_all(indx);
    spiral_left_match_all_perm1 = spiral_left_match_all_perm(indx);
    spiral_right_match_all_perm1 = spiral_right_match_all_perm(indx);
    [edges,ratio_all_mean,ratio_all_sem,ratio_all] = ratio_amp(T2,amp_folder,spiral_left_match_all1,spiral_right_match_all1);
    [edges,ratio_all_mean_perm,ratio_all_sem_perm,ratio_all_perm] = ratio_amp(T2,amp_folder,spiral_left_match_all_perm1,spiral_right_match_all_perm1);    
    [pp(:,count1),p] = anova_test2(edges1,ratio_all,ratio_all_perm);

    errorbar(edges1(1:9),ratio_all_mean(1:9),ratio_all_sem(1:9),color2{kk});
    hold on;
    errorbar(edges1(1:9),ratio_all_mean_perm(1:9),ratio_all_sem_perm(1:9),'k');
    max_ratio = ratio_all_mean+ratio_all_sem;
    
    for i = 1:numel(edges1)
        if p(i)<=0.05 & edges1(i)<=0.025
            text(edges1(1,i),max_ratio(i)+0.02,'*','fontsize',12);
        end
    end  
    hold on;
    yline(0.5,'--k');
    text(0.001,0.85,['amp p = ' num2str(pp(1,count1))],'fontsize',6);
    text(0.001,0.8,['permute p = ' num2str(pp(2,count1))],'fontsize',6);
    xlim([0,0.025]);
    ylim([0.3,0.9]);
    yticks([0.3:0.1:0.9]);
    yticklabels({'0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
    count1 = count1+1;
end
print(h1, 'sprial_matching_ratio_vs_amp_area2', '-dpdf', '-bestfit', '-painters');
%% spirals for each session
count1 = 1;
pixSize = 3.45/1000/0.6*3;  % mm / pix
radius = 40:10:100;
radius2 = radius*pixSize;
p_radius = {};
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    subplot(1,3,count1);
    clear indx T2
    indx = contains(T1.Area,area(1,kk));
    T2 = T1(indx,:);
    spiral_left_match_all1 = spiral_left_match_all(indx);
    spiral_right_match_all1 = spiral_right_match_all(indx);
    spiral_left_match_all_perm1 = spiral_left_match_all_perm(indx);
    spiral_right_match_all_perm1 = spiral_right_match_all_perm(indx);
    [mean_ratio,std_ratio,sem_ratio,ratio_all1] = ratio_radius(T2,spiral_left_match_all1,spiral_right_match_all1);
    [mean_ratio_perm,std_ratio_perm,sem_ratio_perm,ratio_all1_perm] = ratio_radius(T2,spiral_left_match_all_perm1,spiral_right_match_all_perm1);
    
    [pp1(:,count1),p1] = anova_test2(radius,ratio_all1',ratio_all1_perm');
    p_radius{count1} = p1;
    errorbar(radius2,mean_ratio,sem_ratio,color2{kk});
    hold on;
    errorbar(radius2,mean_ratio_perm,sem_ratio_perm,'k');
    hold on;
    yline(0.5,'--k');
    
    max_ratio = mean_ratio+sem_ratio;
    for i = 1:numel(radius)
        if p1(i)<=0.05
            text(radius2(1,i),max_ratio(i)+0.02,'*','fontsize',12);
        end
    end  
    text(0.6,0.68,['radius p = ' num2str(pp1(1,count1))],'fontsize',6);
    text(0.6,0.65,['permute p = ' num2str(pp1(2,count1))],'fontsize',6);
    
    ylim([0.3,0.9]);
    yticks([0.3:0.1:0.9])
    yticklabels({'0.3','0.4','0.5','0.6','0.7','0.8','0.9'})
    count1 = count1+1;
end
%%
print(h1, 'sprial_matching_ratio_neighbor_vs_spiral_size_area', '-dpdf', '-bestfit', '-painters');