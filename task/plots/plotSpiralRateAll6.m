function h5e = plotSpiralRateAll6(data_folder,save_folder)
%%
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
spiral_count_pre_all = [];
spiral_count_post_all = [];
radius = 100;
for i = 1:4
    %%
    clear spiral_count_pre spiral_count_post
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder1,'sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    load(fullfile(data_folder2,[mn '_spirals_task_sort.mat']),'spiral_all');
    load(fullfile(data_folder1,'task_outcome',[mn '_task_outcome']));
    load(fullfile(data_folder1,'task_outcome',[mn '_task_trial_ID']));
    Sessions = unique(trial_all.session);
    %%
    for m = 1:numel(Sessions)
        %%
        trial_index = (trial_all.session == m);
        spiral_session = spiral_all(trial_index,:);
        T_session = T_all(trial_index,:);
        trial_session = trial_all(trial_index,:);
        [spiral_count_pre(m,:), spiral_count_post(m,:)] = getSpiralsCountsBySession2(T_session,spiral_session,radius);
    end
    %%
    spiral_count_pre_all{i,1} = spiral_count_pre;
    spiral_count_post_all{i,1} = spiral_count_post;
end
%%
spiral_count_pre_all1 = cat(1,spiral_count_pre_all{:});
spiral_count_post_all1 = cat(1,spiral_count_post_all{:});
sessionN = cellfun(@(x) size(x,1),spiral_count_pre_all);
%%
ids = [];
sessions_all  = [];
for id = 1:4
    ids = [ids;ones(sessionN(id),1)*id];
    sessions_all = [sessions_all;[1:sessionN(id)]'];
end
ids = repmat(ids,[1,3]);
sessions_all = repmat(sessions_all,[1,3]);
ids = cat(2,ids,ids);
sessions_all = cat(2,sessions_all,sessions_all);
spiral_count = cat(2,spiral_count_pre_all1,spiral_count_post_all1);
labela = ones(35,3);
labels = cat(2,labela,labela*2);
outcome = [1,2,3];
outcome = repmat(outcome,35,1);
outcome = cat(2,outcome,outcome);
%%
Spiral = spiral_count(:);
Session = sessions_all(:);
TimePoint = labels(:);
Subject = ids(:);
Outcome = outcome(:);
Ta = table(Subject,Session,TimePoint,Outcome,Spiral); 
Ta{isnan(Ta{:,5}),5} = 0;
%%
T = readtable('spirals_pre_post2.csv');
T = T(T.TimePoint == 2,:);
T = T(:,[1,2,4,5]);
data.Subject = T.Subject;
data.TimePoint = T.Outcome;
data.Value = T.Value;
data.Sessions = T.Session;
results = permutation_test_paired_timepoints(data);
%%
figure;
cols = colororder;
subplot(1,2,1);
hold on;
for sid = 1:4
    plot(Ta.Spiral(Ta.Subject == sid & Ta.Outcome == 1 & Ta.TimePoint ==2), ...
        Ta.Spiral(Ta.Subject == sid & Ta.Outcome == 2 & Ta.TimePoint ==2), ...
        'o', 'MarkerEdgeColor', cols(sid,:));
    
end
hold on;
plot([0,6]',[0,6]','k--');
axis square;
xlim([0 6]); ylim([0 6]);
xlabel('correct, post');
ylabel('incorrect, post');

subplot(1,2,2);
hold on; 
for sid = 1:4 
    corrPre = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == 1 & Ta.TimePoint == 1);
    corrPost = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == 1 & Ta.TimePoint == 2);
    incorrPre = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == 2 & Ta.TimePoint == 1);
    incorrPost = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == 2 & Ta.TimePoint == 2);
    
    plot(corrPost-corrPre, ...
        incorrPost-incorrPre, ...
        'o', 'MarkerEdgeColor', cols(sid,:));
    
end
hold on;
plot([-3,4]',[-3,4]','k--');
axis square
xlim([-3 4]); ylim([-3 4]);
xlabel('correct, post-pre');
ylabel('incorrect, post-pre');
%%
corrPre = Ta(Ta.outcome_all== "1" & Ta.label_all=="1",:);
corrPost = Ta(Ta.outcome_all== "1" & Ta.label_all=="2",:);
incorrPre = Ta(Ta.outcome_all== "2" & Ta.label_all=="1",:);
incorrPost = Ta(Ta.outcome_all== "2" & Ta.label_all=="2",:);
v = (corrPost{:,5}-corrPre{:,5})-(incorrPost{:,5}-incorrPre{:,5});
tb = corrPre;
tb{:,5} = v;
%%
Ta2 = Ta(Ta.outcome_all == "1"| Ta.outcome_all == "2",:);
writetable(Ta2,'spirals_pre_post2.csv');
%% comapre pre and post in correct trials
for jj = 1:3
    label = string(jj);
    Ta1 = Ta(Ta.outcome_all == label,:);
    lme = fitlme(Ta1, 'spirals ~ label_all +(label_all|id_all)');
    pp(jj,1) = lme.Coefficients{2,6}; 
end    
%% compare correct and incorrect trials in post response
Ta2 = Ta(Ta.label_all == "2",:);
lme2 = fitlme(Ta2, 'spirals ~ outcome_all +(outcome_all|id_all)')
%% compare correct and incorrect trials in pre response
Ta2 = Ta(Ta.label_all == "1",:);
lme3 = fitlme(Ta2, 'spirals ~ outcome_all +(outcome_all|id_all)')
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
