function [T_all] = getTrialResult(T1,data_folder)
T_all = [];
for kk = 1:size(T1,1)
    mn = char(T1.MouseID{kk});
    tda = T1.date(kk);  
    en = T1.folder(kk); 
    block_en = T1.block_folder(kk); 
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load data and block 
    fname = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'task','task_svd',fname);
    expDir = dir(fullfile(session_root,'*_Block.mat'));
    load(fullfile(expDir.folder,expDir.name));
    %% get behavior results
    [T,~] = parse_action_table(block);   
    %%
    ntrial = numel(block.events.endTrialValues);
    norepeatValues = block.events.repeatNumValues(1:ntrial);                 % sometimes, last trial don't have a response
    norepeat_indx = (norepeatValues==1);
    wheel_onset = getWheelOnsetTime(session_root,block);
    wheel_onset = wheel_onset(norepeat_indx);
    %%
    T.wheel_onset = wheel_onset;
    T_all = [T_all;T];
end