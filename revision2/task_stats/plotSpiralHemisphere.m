function h3 = plotSpiralHemisphere(data_folder,save_folder)
%%
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
%%
spirals_all = [];
for i = 1:4
    %%
    clear spiral_count_pre spiral_count_post spiral_count_sum
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder1,'sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    load(fullfile(data_folder2,[mn '_spirals_task_sort.mat']),'spiral_all');
    load(fullfile(data_folder1,'task_outcome',[mn '_task_outcome']));
    load(fullfile(data_folder1,'task_outcome',[mn '_task_trial_ID']));
    Sessions = unique(trial_all.session);
    %%
    spiral = [];
    for m = 1:numel(Sessions)
        trial_index = (trial_all.session == m);
        spiral_session = spiral_all(trial_index,:);
        T_session = T_all(trial_index,:);
        trial_session = trial_all(trial_index,:);
        %%
        for kk = 1:3
            label = labels{kk};
            spiral(m,kk,:,:) = getspiraslHemisphere(T_session,spiral_session,label);
        end    
    end
    spirals_all{i,1} = spiral;
end
%% correct trials
clear spiral_count_sum_contra  spiral_count_sum_ipsi
for i = 1:4
    for k = 1:2
        spirals_temp = spirals_all{i};
        spirals_temp2 = squeeze(spirals_temp(:,1,:,k));
        spirals_mean1 = mean(spirals_temp2,1);
        if k == 1
            spiral_count_sum_contra(i,:) = spirals_mean1;
        else
            spiral_count_sum_ipsi(i,:) = spirals_mean1;
        end
    end
end

t1 = -2:1/35:2;
indx = find(t1>=-1 & t1<=1);
t1 = t1(indx);
spiral_count_sum_contra = spiral_count_sum_contra(:,indx);
spiral_count_sum_ipsi = spiral_count_sum_ipsi(:,indx);
ylim2 = 10;
mouseN = 4;
h2 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
subplot(1,2,1);
spiral_mean_contra = mean(spiral_count_sum_contra,1);
spiral_sem_contra = std(spiral_count_sum_contra,1)./sqrt(mouseN);
shadedErrorBar(t1, spiral_mean_contra, spiral_sem_contra, 'lineprops', 'g');
xline(t1(36),'k--');
xline(t1(39),'k--');
xline(t1(46),'k--');
ylim([2,ylim2]);
ylabel('Rotating waves/s');

subplot(1,2,2);
spiral_mean_ipsi = mean(spiral_count_sum_ipsi,1);
spiral_sem_ipsi = std(spiral_count_sum_ipsi,1)./sqrt(mouseN);
shadedErrorBar(t1, spiral_mean_ipsi, spiral_sem_ipsi, 'lineprops', 'g');
xline(t1(36),'k--');
xline(t1(39),'k--');
xline(t1(46),'k--');
ylim([2,ylim2]);
ylabel('Rotating waves/s');
print(h2,'spirals_left_right_hemisphere.pdf','-dpdf', '-bestfit', '-painters');
%% stats
frames_pre = 63:67; frames_post = 77:81;
for i = 1:4
    spirals_temp = spirals_all{i};
    spiral_count_pre = squeeze(sum(spirals_temp(:,:,frames_pre,2),3))/numel(frames_pre);
    spiral_count_post = squeeze(sum(spirals_temp(:,:,frames_post,2),3))/numel(frames_pre);
    %%
    spiral_count_pre_all{i,1} = spiral_count_pre;
    spiral_count_post_all{i,1} = spiral_count_post;
end

spiral_count_pre_all1 = cat(1,spiral_count_pre_all{:});
spiral_count_post_all1 = cat(1,spiral_count_post_all{:});
sessionN = cellfun(@(x) size(x,1),spiral_count_pre_all);

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
Spiral = spiral_count(:);
Session = sessions_all(:);
TimePoint = labels(:);
Subject = ids(:);
Outcome = outcome(:);
Ta = table(Subject,Session,TimePoint,Outcome,Spiral); 
Ta{isnan(Ta{:,5}),5} = 0;
%% compare pre and post
for k = 1:3
    T = Ta(Ta.Outcome == k,:);
    data.Subject = T.Subject;
    data.TimePoint = T.TimePoint;
    data.Value = T.Spiral;
    data.Sessions = T.Session;
    results = permutation_test_paired_timepoints(data,'makePlots',false,'showProgress',false);
    pvalues(k,1)= results.p_value_two_tailed;
end

% compare evoked spirals amongst different outcome types
% [correct vs incorrect ; correct vs miss; incorrect vs miss]
exclude_outcome = [3,2,1];
for k = 1:3
    T = Ta(not(Ta.Outcome == exclude_outcome(k)),:);
    if k == 2
        T.Outcome(T.Outcome==3) = 2;
    elseif k == 3
        T.Outcome(T.Outcome==3) = 1;
    end
    diff = T{T.TimePoint == 2,5}-T{T.TimePoint == 1,5}; % post-pre
    T = T(T.TimePoint == 2,:); %use labels at post response
    T = T(:,[1,2,4,5]);
    data.Subject = T.Subject; % data structure
    data.TimePoint = T.Outcome; % rename outcome to timepoints for function variable convinience
    data.Value = diff; % use diff as value
    data.Session = T.Session;
    %%
    results = permutation_test_paired_timepoints(data,'makePlots',false,'showProgress',false);
    pvalues2(k,1)= results.p_value_two_tailed;
end
%%
outcome2 = [1,2;1,3;2,3];
symbols = {"o","+","x","square"};
axis_label2 = {"Correct,pre","Correct,post";"Incorrect,pre","Incorrect,post";"Miss,pre","Miss,post"};
axis_label = {"Incorrect","Correct";"Miss","Correct";"Miss","Incorrect"};
h2 = figure('Renderer', 'painters', 'Position', [50 50 1000 600]);
cols = cbrewer2('qual','Set1',4);
for k = 1:3
    subplot(2,3,k);
    hold on;
    for sid = 1:4
        pre = Ta.Spiral(Ta.Subject == sid & Ta.Outcome == k & Ta.TimePoint == 1);
        post = Ta.Spiral(Ta.Subject == sid & Ta.Outcome == k & Ta.TimePoint == 2);
        p = plot(pre,post,'Marker',symbols{sid}, 'MarkerEdgeColor', 'k','MarkerSize', 4);
        p.LineStyle = 'none';
        %%
        pre_mean = mean(pre); pre_sem = std(pre)./sqrt(numel(pre));
        post_mean = mean(post); post_sem = std(post)./sqrt(numel(post));
        p2 = plot(pre_mean,post_mean,'Marker',symbols{sid}, 'MarkerEdgeColor', cols(sid,:),...
            'MarkerSize', 4,'LineWidth',1);
        hold on
        eb(1) = errorbar(pre_mean,post_mean,pre_sem, 'horizontal', 'LineStyle', 'none');
        eb(2) = errorbar(pre_mean,post_mean,post_sem, 'vertical', 'LineStyle', 'none');
        set(eb, 'color', cols(sid,:), 'LineWidth', 0.5)
        %%
        pre_mean_all(sid,k) = pre_mean;
        post_mean_all(sid,k) = post_mean;
    end
    text(10,30,['p = ' num2str(round(pvalues(k)*10000)/10000)]);
    hold on;
    plot([0,15]',[0,15]','k--');
    xlim([0 15]); ylim([0 15]);
    axis square;
    xlabel(axis_label2(k,1));
    ylabel(axis_label2(k,2));
    %%
    subplot(2,3,k+3);
    hold on; 
    for sid = 1:4 
        typea_pre = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == outcome2(k,1) & Ta.TimePoint == 1);
        typea_post = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == outcome2(k,1) & Ta.TimePoint == 2);
        typeb_pre = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == outcome2(k,2) & Ta.TimePoint == 1);
        typeb_post = Ta.Spiral(Ta.Subject ==sid & Ta.Outcome == outcome2(k,2) & Ta.TimePoint == 2);
        typea_diff = typea_post-typea_pre;
        typeb_diff = typeb_post-typeb_pre;
        p = plot(typeb_diff, typea_diff, ...
            'Marker',symbols{sid}, 'MarkerEdgeColor','k','MarkerSize', 4);
        p.LineStyle = 'none';
        %%
        typea_mean = mean(typea_diff); typea_sem = std(typea_diff)./sqrt(numel(typea_diff));
        typeb_mean = mean(typeb_diff); typeb_sem = std(typeb_diff)./sqrt(numel(typeb_diff));
        p2 = plot(typeb_mean,typea_mean,'Marker',symbols{sid}, 'MarkerEdgeColor', cols(sid,:),...
            'MarkerSize', 4,'LineWidth',1);
        hold on
        eb(1) = errorbar(typeb_mean,typea_mean,typeb_sem, 'horizontal', 'LineStyle', 'none');
        eb(2) = errorbar(typeb_mean,typea_mean,typea_sem, 'vertical', 'LineStyle', 'none');
        set(eb, 'color', cols(sid,:), 'LineWidth', 0.5)
        %%
        typea_diff_all(sid,k) = typea_mean;
        typeb_diff_all(sid,k) = typeb_mean;
    end
    %%
    text(-8,16,['p = ' num2str(round(pvalues2(k)*10000)/10000)]);
    hold on;
    plot([-16,16]',[-16,16]','k--');
    axis square
    xlim([-16 16]); ylim([-16 16]);
    axis_lim = [-16:4:16];
    xticks(axis_lim); yticks(axis_lim);
    xticklabels(axis_lim); yticklabels(axis_lim);
    xlabel(axis_label(k,1));
    ylabel(axis_label(k,2));
end