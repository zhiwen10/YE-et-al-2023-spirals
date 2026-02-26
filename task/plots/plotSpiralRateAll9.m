function h5e = plotSpiralRateAll9(data_folder,save_folder)
%%
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
spiral_count_sum_all = [];
spiral_count_pre_all = [];
spiral_count_post_all = [];
% radii = [40:10:100];
radii = [40:10:100];
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
    for m = 1:numel(Sessions)
        clear spiral_count_sum_temp
        %%
        trial_index = (trial_all.session == m);
        spiral_session = spiral_all(trial_index,:);
        T_session = T_all(trial_index,:);
        trial_session = trial_all(trial_index,:);
        for r = 1:numel(radii)
            radius = radii(r);
            spiral_count_sum_temp = getSpiralsCountsBySession3(T_session,spiral_session,radius);
            spiral_count_sum(m,r,:,:) = spiral_count_sum_temp;
            spiral_count_pre(m,r,:) = squeeze(sum(spiral_count_sum_temp(:,64:66),2))/3;
            spiral_count_post(m,r,:) = squeeze(sum(spiral_count_sum_temp(:,78:80),2))/3;
        end
    end
    %%
    spiral_count_sum_all{i,1} = spiral_count_sum;
    spiral_count_pre_all{i,1} = spiral_count_pre;
    spiral_count_post_all{i,1} = spiral_count_post;
end
%%
spiral_count_pre_all1 = cat(1,spiral_count_pre_all{:});
spiral_count_post_all1 = cat(1,spiral_count_post_all{:});
sessionN = cellfun(@(x) size(x,1),spiral_count_pre_all);
%%
sessionN2 = cumsum(sessionN);
sessionN2 = [0;sessionN2];
%%
for id = 1:4
    %%
    spiral_pre_temp = spiral_count_pre_all1(sessionN2(id)+1:sessionN2(id+1),:,:);
    spiral_post_temp = spiral_count_post_all1(sessionN2(id)+1:sessionN2(id+1),:,:);
    spiral_pre_mean(id,:,:) = mean(spiral_pre_temp);
    spiral_post_mean(id,:,:) = mean(spiral_post_temp);
end
%%
spiral_pre_mean2 = squeeze(mean(spiral_pre_mean,1));
spiral_pre_sem2 = squeeze(std(spiral_pre_mean,1))./sqrt(4);
spiral_post_mean2 = squeeze(mean(spiral_post_mean,1));
spiral_post_sem2 = squeeze(std(spiral_post_mean,1))./sqrt(4);
%%
for r = 1:7
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
    spiral_count_pre_all1a = squeeze(spiral_count_pre_all1(:,r,:));
    spiral_count_post_all1a = squeeze(spiral_count_post_all1(:,r,:));
    spiral_count = cat(2,spiral_count_pre_all1a,spiral_count_post_all1a);
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
    % compare pre and post
    for k = 1:3
        T = Ta(Ta.Outcome == k,:);
        data.Subject = T.Subject;
        data.TimePoint = T.TimePoint;
        data.Value = T.Spiral;
        data.Sessions = T.Session;
        results = permutation_test_paired_timepoints(data,'makePlots',false,'showProgress',false);
        pvalues(k,r)= results.p_value_two_tailed;
    end
end
%%
for k = 1:3
    [corrected_p(k,:), hha(k,:)]=bonf_holm(pvalues(k,:),0.05);
end
%%
t1 = -2:1/35:2;
indx = find(t1>=-1 & t1<=1);
t2 = t1(indx);
cols = cbrewer2('qual','Set1',4);
h2 = figure('Renderer', 'painters', 'Position', [50 50 1000 600]);
for k = 1:3
    for r = 1:7
        for id = 1:4
            trace_temp = spiral_count_sum_all{id};
            trace_temp2 = squeeze(trace_temp(:,:,k,:));
            trace_temp2a = squeeze(trace_temp2(:,r,:));
            trace_mean(id,:) = mean(trace_temp2a,1,'omitnan');  
        end
        %%
        subplot(3,7,r+(k-1)*7);                      
        trace_mean2 = mean(trace_mean,1,'omitnan');
        trace_sem2 = std(trace_mean,1,'omitnan')./sqrt(size(trace_mean,1));
        shadedErrorBar(t2, trace_mean2(indx),trace_sem2(indx), 'lineprops', 'g');
        xline(0,'k--');
        xline(0.25,'k--');
        ylim([0,5]);
        xticks([-1:0.5:1]);
    end
end
%%
print(h2,'spirals_all_radius_timeline.pdf','-dpdf', '-bestfit', '-painters');