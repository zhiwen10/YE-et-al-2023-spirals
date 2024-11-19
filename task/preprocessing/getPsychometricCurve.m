function getPsychometricCurve(data_folder,save_folder)
%%
fname_all = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
rt_median_all_mean = nan(11,4); rt_median_all_sem = nan(11,4);
left_choice_all_mean = nan(11,4); left_choice_all_sem = nan(11,4);
right_choice_all_mean = nan(11,4); right_choice_all_sem = nan(11,4);
no_go_all_mean = nan(11,4); no_go_all_sem = nan(11,4);
T_trial_counts_all = {};
%%
for m = 1:4
    %%
    clearvars -except rt_median_all_mean rt_median_all_sem left_choice_all_mean...
        left_choice_all_sem right_choice_all_mean  right_choice_all_sem ...
        no_go_all_mean no_go_all_sem m fname_all T_trial_counts_all data_folder save_folder
    T_session = readtable(fullfile(data_folder,'task','sessions', [fname_all{m} '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    sessions = size(T1,1);
    rt_median_all = nan(11,sessions);
    left_choice_all = nan(11,sessions);
    right_choice_all = nan(11,sessions);
    no_go_all = nan(11,sessions);
    %%
    for kk = 1:sessions
        %%
        clearvars -except T_session T1 sessions T_trial_counts_all kk m...
            rt_median_all left_choice_all right_choice_all no_go_all...
            rt_median_all_mean rt_median_all_sem left_choice_all_mean...
            left_choice_all_sem right_choice_all_mean  right_choice_all_sem ...
            no_go_all_mean no_go_all_sem fname_all data_folder save_folder
        mn = char(T1.MouseID{kk});
        tda = T1.date(kk);  
        en = T1.folder(kk); 
        block_en = T1.block_folder(kk); 
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        %%   
        % load block
        fname = [mn '_' tdb '_' num2str(en)];
        session_root = fullfile(data_folder,'task','task_svd',fname);
        load(fullfile(session_root,[td '_' num2str(block_en) '_' mn '_Block.mat']));
        ntrial = numel(block.events.endTrialValues);
        norepeatValues = block.events.repeatNumValues(1:ntrial);                 % sometimes, last trial don't have response
        response = block.events.responseValues(1:ntrial);
        norepeat_indx = (norepeatValues==1);
        % get response no repeat
        response = response(norepeat_indx)';
        % get feedback no repeat
        feedback = block.events.feedbackValues(norepeat_indx)';
        % get stim no repeat
        left_contrast = block.events.contrastLeftValues;
        right_contrast = block.events.contrastRightValues;
        left_contrast = left_contrast(norepeat_indx)';
        right_contrast = right_contrast(norepeat_indx)';
        % reactiom time
        reaction_time = block.events.responseTimes(1:ntrial)-block.events.stimulusOnTimes(1:ntrial);
        reaction_time = reaction_time(norepeat_indx)';
        %%
        wheel_onset = getWheelOnsetTime(session_root,block);
        wheel_onset = wheel_onset(norepeat_indx);
        %% get table
        T = table(left_contrast,right_contrast,response,feedback);
        % miss
        miss = zeros(size(T,1),1);
        miss(logical(T.right_contrast-T.left_contrast) & T.response == 0,1)= 1;
        T.miss = miss;
        % correct
        correct = zeros(size(T,1),1);
        correct((logical(T.right_contrast-T.left_contrast)& T.feedback == 1),1)= 1;
        T.correct = correct;
        % incorrect
        incorrect = zeros(size(T,1),1);
        % incorrect((logical(contrastRight-contrastLeft)& not(miss) & (choice~=action)),1)= 1;
        incorrect((logical(T.right_contrast-T.left_contrast)& logical(T.response) & T.feedback==0),1)= 1;
        T.incorrect = incorrect;
        % reject
        reject = zeros(size(T,1),1);
        reject((T.right_contrast == 0 & T.left_contrast == 0 & T.response==0),1)= 1;
        T.reject = reject;
        % falarm
        falarmL = zeros(size(T,1),1);
        falarmL((T.right_contrast == 0 & T.left_contrast == 0& T.response == 1),1) = 1;
        T.falarmL = falarmL;

        falarmR = zeros(size(T,1),1);
        falarmR((T.right_contrast == 0 & T.left_contrast == 0& T.response == -1),1) = 1;
        T.falarmR = falarmR;
        % check if all trials were conted

        sum_T = sum(T{:,5:10},2);
        match = all(sum_T);
        T_label = T(:,5:10);
        label1 = {};
        list = {'miss','correct','incorrect','reject','falarmL','falarmR'};
        for i = 1:size(T_label,1)
            binary_temp = T_label{i,:};
            a = list{find(binary_temp)};
            label1(i,1) = {a};
        end
        T.label = label1;
        T.reactionTime = reaction_time;
        T.wheel_onset = wheel_onset;
        %%
        T.left_contrast(T.left_contrast==0.006) = 0.06;
        T.right_contrast(T.right_contrast==0.006) = 0.06;
        %%
        contrast = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1]'; 
        stim = [1, 0; 0.5, 0; 0.25, 0;0.125, 0; 0.06,0; 0, 0; 0,0.06; 0,0.125; 0, 0.25; 0, 0.5; 0, 1];
        ratio = [];
        trial_counts = [];
        reaction_time_sort = {};
        T_sort = {};
         for i = 1:11
             clear T_temp
             T_temp = T(T.left_contrast==stim(i,1) & T.right_contrast==stim(i,2),:);
             ratio(i,:) = sum(T_temp{:,5:10},1)./size(T_temp,1);
             trial_counts(i,:) = sum(T_temp{:,5:10},1);  
             sum_counts(i,:) = size(T_temp,1);
             T_sort{i} = T_temp;
             reaction_time_sort{i} = T_temp(ismember(T_temp.label,{'correct','falarmL','falarmR'}),:).wheel_onset;
             % reaction_time_sort{i} = T_temp(ismember(T_temp.label,{'correct','falarmL','falarmR'}),:).reactionTime;
         end
        %%
        rt_median = cellfun(@median, reaction_time_sort)';
        %%
        ratio = [contrast,ratio];
        T_ratio = array2table(ratio);
        T_ratio = renamevars(T_ratio,["ratio1","ratio2","ratio3","ratio4","ratio5","ratio6","ratio7"], ...
                         ["contrast","miss","correct","incorrect","reject","falarmL","falarmR"]);
        trial_counts = [contrast, trial_counts, sum_counts];
        T_trial_counts = array2table(trial_counts);
        T_trial_counts = renamevars(T_trial_counts,["trial_counts1","trial_counts2","trial_counts3",...
            "trial_counts4","trial_counts5","trial_counts6","trial_counts7","trial_counts8"], ...
            ["contrast","miss","correct","incorrect","reject","falarmL","falarmR","total"]);           
        left_choice = [T_ratio{1:5,4};T_ratio{6,6};T_ratio{7:11,3}];
        right_choice = [T_ratio{1:5,3};T_ratio{6,7};T_ratio{7:11,4}];
        no_go = [T_ratio{1:5,2};T_ratio{6,5};T_ratio{7:11,2}];
        %%
        T_trial_counts_all{kk,m} = T_trial_counts;
        rt_median_all(:,kk) = rt_median;
        left_choice_all(:,kk) = left_choice;
        right_choice_all(:,kk) = right_choice;
        no_go_all(:,kk) = no_go;
    end
    %%
    rt_median_all_mean(:,m) = mean(rt_median_all,2,'omitnan');
    rt_median_all_sem(:,m) = std(rt_median_all,[],2,'omitnan')./sqrt(size(rt_median_all,2));
    left_choice_all_mean(:,m) = mean(left_choice_all,2,'omitnan');
    left_choice_all_sem(:,m) = std(left_choice_all,[],2,'omitnan')./sqrt(size(left_choice_all,2));
    right_choice_all_mean(:,m) = mean(right_choice_all,2,'omitnan');
    right_choice_all_sem(:,m) = std(right_choice_all,[],2,'omitnan')./sqrt(size(right_choice_all,2));
    no_go_all_mean(:,m) = mean(no_go_all,2,'omitnan');
    no_go_all_sem(:,m) = std(no_go_all,[],2,'omitnan')./sqrt(size(no_go_all,2));
end
%%
save(fullfile(save_folder,'Task_performance_across_sessions'),...
    'T_trial_counts_all','rt_median_all_mean','rt_median_all_sem',...
    'left_choice_all_mean','left_choice_all_sem','right_choice_all_mean','right_choice_all_sem',...
    'no_go_all_mean','no_go_all_sem');