function h3 = plotSpiralDirection2(data_folder,save_folder)
%%
data_folder = 'E:\spiral_data_share\data';  
% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
scale = 1;
[row,col] = find(BW);
brain_index = [col,row];
%% right SSp index
clear areaPath
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
Utransformed = projectedAtlas1;
scale = 1;

hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp_right = [col,row];

hemi = 'left';
[indexSSp2,UselectedSSp2] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexSSp2);
indexSSp_left = [col2,row2];
%%
% plot task trials with spiral raster
data_folder1 = fullfile(data_folder,'task');
data_folder2 = fullfile(data_folder,'task','spirals');
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct","incorrect","miss"};
%%
radius = 100;
ratio_right = [];
ratio_left = [];
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
    ratio_right = [];
    ratio_left = [];
    for m = 1:numel(Sessions)
        trial_index = (trial_all.session == m);
        spiral_session = spiral_all(trial_index,:);
        T_session = T_all(trial_index,:);
        trial_session = trial_all(trial_index,:);
        %%
        for kk = 1:3
            label = labels{kk};
            [ratio_right(m,kk,:),ratio_left(m,kk,:)] = getspiraslHemisphere3...
                (T_session,spiral_session,label,radius,indexSSp_right,indexSSp_left);
        end    
    end
    ratio_contra_all{i,1} = ratio_right;
    ratio_ipsi_all{i,1} = ratio_left;
end
%%
for i = 1:4
    ratio_contra = ratio_contra_all{i};
    ratio_contra2(i,:) = squeeze(mean(ratio_contra(:,1,:),1));
    ratio_ipsi = ratio_ipsi_all{i};
    ratio_ipsi2(i,:) = squeeze(mean(ratio_ipsi(:,1,:),1));
end
%% permutation stats
ratio_contra_all2 = cat(1,ratio_contra_all{:});
ratio_contra_all3 = cat(2,ratio_contra_all2(:,:,1),ratio_contra_all2(:,:,2));
ratio_ipsi_all2 = cat(1,ratio_ipsi_all{:});
ratio_ipsi_all3 = cat(2,ratio_ipsi_all2(:,:,1),ratio_ipsi_all2(:,:,2));
%%
sessionN = cellfun(@(x) size(x,1),ratio_contra_all);
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
labela = ones(35,3);
labels = cat(2,labela,labela*2);
outcome = [1,2,3];
outcome = repmat(outcome,35,1);
outcome = cat(2,outcome,outcome);
%%
Spiral = ratio_contra_all3(:);
% Spiral = ratio_ipsi_all3(:);
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
    % diff = T{T.TimePoint == 2,5}-T{T.TimePoint == 1,5}; % post-pre
    % diff = T{T.TimePoint == 2,5}; % post
    diff = T{T.TimePoint == 1,5}; % pre
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
    text(0,1,['p = ' num2str(round(pvalues(k)*10000)/10000)]);
    hold on;
    plot([0,1]',[0,1]','k--');
    xlim([0 1]); ylim([0 1]);
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
%         typea_diff = typea_post-typea_pre;
%         typeb_diff = typeb_post-typeb_pre;
%         typea_diff = typea_post;
%         typeb_diff = typeb_post;
        typea_diff = typea_pre;
        typeb_diff = typeb_pre;
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
    text(0,1,['p = ' num2str(round(pvalues2(k)*10000)/10000)]);
    hold on;
    plot([0,1]',[0,1]','k--');
    axis square
    xlim([0 1]); ylim([0 1]);
    axis_lim = [0:0.5:1];
    xticks(axis_lim); yticks(axis_lim);
    xticklabels(axis_lim); yticklabels(axis_lim);
    xlabel(axis_label(k,1));
    ylabel(axis_label(k,2));
end
%%
print(h2,'CW_contra_ratio_pre_large.pdf','-dpdf', '-bestfit', '-painters');
%%
for k = 1:3
    [hh1(k,1),pp1(k,1)] = ttest(pre_mean_all(:,k),post_mean_all(:,k));
    [hh2(k,1),pp2(k,1)] = ttest(typea_diff_all(:,k),typeb_diff_all(:,k));
end
%% scatter plot
outcome2 = [1,2;1,3;2,3];
symbols = {"o","+","x","square"};
axis_label2 = {"Correct,pre","Correct,post";"Incorrect,pre","Incorrect,post";"Miss,pre","Miss,post"};
axis_label = {"Incorrect","Correct";"Miss","Correct";"Miss","Incorrect"};
h2 = figure('Renderer', 'painters', 'Position', [50 50 500 200]);
cols = cbrewer2('qual','Set1',4);
for k = 1:3
    subplot(1,3,k);
    hold on;
    for sid = 1:4
        pre = Ta.Spiral(Ta.Subject == sid & Ta.Outcome == k & Ta.TimePoint == 1);
        post = Ta.Spiral(Ta.Subject == sid & Ta.Outcome == k & Ta.TimePoint == 2);
        scatter(ones(size(pre)),pre,16,'k');
        hold on;
        scatter(ones(size(post))*2,post,16,'k');
        hold on;
        plot([ones(size(pre)),ones(size(post))*2]',[pre,post]','Color',[0.5,0.5,0.5]);
        %%
        pre_mean = mean(pre); pre_sem = std(pre)./sqrt(numel(pre));
        post_mean = mean(post); post_sem = std(post)./sqrt(numel(post));
        scatter(ones(size(pre_mean)),pre_mean,16,'r','filled');
%         hold on;
%         eh1 = errorbar(1,pre_mean,pre_sem, 'vertical', 'Color',cols(sid,:),'LineWidth',2);
%         eh1.CapSize = 12;
        hold on;
        scatter(ones(size(post_mean))*2,post_mean,16,'r','filled');
%         hold on;
%         eh2 = errorbar(2,post_mean,post_sem, 'vertical', 'Color',cols(sid,:),'LineWidth',2);
%         eh2.CapSize = 12;
        hold on;
        plot([ones(size(pre_mean)),ones(size(post_mean))*2]',[pre_mean,post_mean]',...
            'Color','k','LineWidth',2);
        %%
        pre_mean_all(sid,k) = pre_mean;
        post_mean_all(sid,k) = post_mean;
    end
    text(1.2,1,['p = ' num2str(round(pvalues(k)*10000)/10000)]);
    hold on;
    % plot([0,1]',[0,1]','k--');
    yline(0.5,'k--');
    ylim([0 1]);
    % axis square;
    xticks([1,2]);
    xticklabels({'Pre','Post'});
    ylabel('CW rotating wave ratio');
end
%%
print(h2,'CW_contra_ratio_post_large_2.pdf','-dpdf', '-bestfit', '-painters');