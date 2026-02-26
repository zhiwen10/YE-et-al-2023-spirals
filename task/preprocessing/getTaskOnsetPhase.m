function getTaskOnsetPhase(data_folder,save_folder,freq)
%%
fname = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
for i = 1:4
    clear wf_all phase_all amp_all contrast_all
    mn = fname{i};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    win = [0,5000];
    trialWin = [-4,4];
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
    %%
    [wf_all,filt_all,phase_all,amp_all, contrast_all] = ...
        getTrialTraceTask3(data_folder,T1,win,trialWin, freq);
    save(fullfile(save_folder,[mn '_task_freq_to' num2str(freq(2)) 'Hz']),...
        'wf_all','filt_all', 'phase_all','amp_all','contrast_all')
end
