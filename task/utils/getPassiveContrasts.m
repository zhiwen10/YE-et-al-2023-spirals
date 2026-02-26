function [contrast_all] = getPassiveContrasts(T1,data_folder,win)
contrast_all = [];
for kk = 1:size(T1,1)
    %%
    clear contrasts
    mn = char(T1.MouseID{kk});
    tda = T1.date(kk);  
    en = T1.folder(kk); 
    block_en = T1.block_folder(kk); 
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    % load data and block 
    fname = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'task','task_svd',fname);
    expDir = dir(fullfile(session_root,'*_Block.mat'));
    load(fullfile(expDir.folder,expDir.name));
    % get photodiode time
    allPD2 = getPhotodiodeTime2(session_root,block,win);
    ntrials = numel(block.events.endTrialValues);
    left_contrast = block.events.contrastLeftValues(1:ntrials);
    right_contrast = block.events.contrastRightValues(1:ntrials);
    contrasts = [left_contrast; right_contrast]';
    %
    contrast_all = [contrast_all;contrasts];
end