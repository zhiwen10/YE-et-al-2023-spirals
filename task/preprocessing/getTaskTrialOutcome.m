function getTaskTrialOutcome(data_folder, save_folder)
%%
fname = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
for i = 1:4
    %%
    clear T_all
    mn = fname{i};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    [T_all] = getTrialResult(T1,data_folder);    
    save(fullfile(save_folder,[mn '_task_outcome']),'T_all');
    [trial_all] = getTrialID(T1,data_folder);
    save(fullfile(save_folder,[mn '_task_trial_ID']),'trial_all');
end
