githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\AP_scripts_cortexlab-master'))
addpath(genpath(fullfile(githubDir, 'wheelAnalysis')))
addpath(genpath(fullfile(githubDir, 'PatternDetection')))
addpath(genpath(fullfile(githubDir, 'spiralDetection')));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\shadedErrorBar'));
%%
area{1} = 'THAL'; area{2} = 'STR'; area{3} = 'CORTEX'; area{4} = 'MB';
%% total sessions
codefolder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow1';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T1  = T(T.meanVar>0&T.maxVar>0 & not(isnan(T.meanAmp)),:);
T1(contains(T1.Area, 'CORTEX'),:) = [];
count1 = 1;
for current_area = [1,2,4]
    indx = contains(T.Area,area{current_area});
    total_count(count1) = sum(indx);
    count1 = count1+1;
end
%% mean and max variances
max_var = {}; mean_var = {};
count1 = 1;
for current_area = [1,2,4]
    indx = find(contains(T1.Area,area{current_area}) & T1.meanVar>=0.1);
    current_T = T1(indx,:);
    list = 1:size(current_T,1);    
    max_var{count1} = current_T.maxVar;
    mean_var{count1} = current_T.meanVar;
    count1 = count1+1;
end
% significance test
mean_var_all = cellfun(@mean,mean_var);
max_var_all = cellfun(@mean,max_var);
std_mean_all = cellfun(@std,mean_var);
std_max_all = cellfun(@std,max_var);
count = cellfun(@numel,max_var);
sem_mean_all = std_mean_all./sqrt(count);
sem_max_all = std_max_all./sqrt(count);
pairs = {[1,2],[1,3],[2,3]};
for i = 1:3
    mean_pair = mean_var(pairs{i});
    max_pair = max_var(pairs{i});
    [h1(i,1),p1(i,1)] = ttest2(mean_pair{1},mean_pair{2});
    [h1(i,2),p1(i,2)] = ttest2(max_pair{1},max_pair{2});
end
Ta1 = array2table(p1,...
    'VariableNames',{'mean_var','max_var'});
Ta1.name = {'thal vs str'; 'thal vs mb'; 'str vs mb'};
Ta1 = [Ta1(:,3),Ta1(:,1:2)];
%%
writetable(Ta1,'mean_max_var_significance_test.csv');
%% neuron counts
data_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\variance_summary\variance_order\var_ordered';
count1 = 1;
numVar_all = {};
for current_area = [1,2,4]
    indx = find(contains(T1.Area,area{current_area})  & T1.meanVar>=0.1);
    current_T = T1(indx,:);
    %%
    numVar = [];
    for kk = 1:size(current_T,1) 
        ops = get_session_info(current_T,kk);
        var_name = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_explained_var.mat'];
        load(fullfile(data_folder,var_name),'mean_var2');
        max_var1 = max(mean_var2);
        indx = find(mean_var2>=0.8*max_var1,1,'first');
        var_current = mean_var2(indx);
        numVar(kk,:) = [numel(mean_var2),indx,var_current];
    end
    numVar_all{count1} = numVar;
    count1 = count1+1;
end
mean_area = cellfun(@mean,numVar_all,'UniformOutput',false);
mean_area = cat(1,mean_area{:});
std_area = cellfun(@std,numVar_all,'UniformOutput',false);
std_area = cat(1,std_area{:});
count_area = cellfun(@(x) size(x,1),numVar_all);
count_area = repmat(count_area',1,3);
sem_area = std_area./sqrt(count_area);
% stats
pairs = {[1,2],[1,3],[2,3]}; % thal vs str; thal vs mb; str vs mb
for i = 1:3
    mean_pair = numVar_all(pairs{i});
    [h(i,1),p(i,1)] = ttest2(mean_pair{1}(:,1),mean_pair{2}(:,1));
    [h(i,2),p(i,2)] = ttest2(mean_pair{1}(:,2),mean_pair{2}(:,2));
    [h(i,3),p(i,3)] = ttest2(mean_pair{1}(:,3),mean_pair{2}(:,3));
end
Ta = array2table(p,'VariableNames',{'total_neurons','Effective_neurons','variances'});
Ta.name = {'thal vs str'; 'thal vs mb'; 'str vs mb'};
Ta = [Ta(:,4),Ta(:,1:3)];
%%
writetable(Ta,'significance_test.csv');
%%
color1 = {'g','r','m'};
count1 = 1;
h1 = figure('Renderer', 'painters', 'Position', [100 100 700 700]);
subplot(2,2,1);
% scatter(T1.maxVar,T1.meanVar,8,'k');
% hold on;
for i = 1:3
    scatter(max_var{i},mean_var{i},8,color1{i},'filled');
    hold on;
    errorbar(max_var_all(i),mean_var_all(i),sem_mean_all(i),sem_mean_all(i),sem_max_all(i),sem_max_all(i),...
    'marker','o','color',color1{i});
    hold on;
end
% axis equal;
xlim([0.2,1.0]); ylim([0,0.8]); 
yticks([0:0.2:0.8]);
yticklabels(string(0:0.2:0.8));
xlabel('Max var explained');
ylabel('Mean var explained');
subplot(2,2,2);
for i = 1:3
    current_numVar = numVar_all{i};
    scatter(current_numVar(:,1),current_numVar(:,3),8,color1{i},'filled');
    hold on;
    errorbar(mean_area(i,1),mean_area(i,3),sem_area(i,3),sem_area(i,3),sem_area(i,1),sem_area(i,1),...
        'marker','o','color',color1{i});
    hold on;
end
% xlim([0,40]);
ylim([0,0.8]);
% xticks([0,10,20,30,40]);
% xticklabels({'0','10','20','30','40'});
yticks([0:0.2:0.8]);
yticklabels(string(0:0.2:0.8));
xlabel('Neuron number');
ylabel('Mean var explained'); 
%
count1 = 1;
for current_area = [1,2,4]
    %%
    subplot(2,6,count1+6);
    indx = find(contains(T1.Area,area{current_area})& T1.meanVar>=0.1);
    current_T = T1(indx,:);
    %%
    for kk = 1:size(current_T,1) 
        %%
        ops = get_session_info(current_T,kk);
        var_name = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_explained_var.mat'];
        load(fullfile(data_folder,var_name),'mean_var2');
        plot(1:numel(mean_var2),mean_var2,'color',color1{count1});
        title(area{current_area});
        hold on;
    end
    xlim([0,200]);
    ylim([0,0.8]);
    xticks([0, 50, 100, 150,200]);
    xticklabels({'0','50','100','150','200'});
    yticks([0:0.2:0.8]);
    yticklabels(string(0:0.2:0.8));
    xlabel('Neuron number');
    ylabel('Mean var explained');
    count1 = count1+1;
end
%
subplot(2,2,4);
hold off;
for i = 1:3
    current_numVar = numVar_all{i};
    scatter(current_numVar(:,2),current_numVar(:,3),8,color1{i},'filled');
    hold on;
    errorbar(mean_area(i,2),mean_area(i,3),sem_area(i,3),sem_area(i,3),sem_area(i,2),sem_area(i,2),...
        'marker','o','color',color1{i});
    hold on;
end
xlim([0,50]);
ylim([0,0.8]);
xticks([0,10,20,30,40,50]);
xticklabels({'0','10','20','30','40','50'});
yticks([0:0.2:0.8]);
yticklabels(string(0:0.2:0.8));
xlabel('Neuron number');
ylabel('80% of mean var explained'); 
%%
print(h1,'variance_summary_new' , '-dpdf', '-bestfit', '-painters');