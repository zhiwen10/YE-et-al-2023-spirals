function getTaskSpirals(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
scale = 1;
[row,col] = find(BW);
brain_index = [col,row];
%%
spiral_folder = fullfile(data_folder,'task','spirals_all');
%%
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
for i = 1:4
    %%
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    sessions = size(T1,1);
    %
    ntrials_all = 0;
    spiral_all = {};
    for kk = 1:size(T1,1)
        clear T T_ratio allPD1 allPD2 stimOn frames_stimOn pwAll spiral_cell spiral_cat
        mn = char(T1.MouseID{kk});
        tda = T1.date(kk);  
        en = T1.folder(kk); 
        block_en = T1.block_folder(kk); 
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        %% load data and block 
        fname = [mn '_' tdb '_' num2str(en)];
        session_root = fullfile(data_folder,'task','task_svd',fname);
        [U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
        dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
        expDir = dir(fullfile(session_root,'*_Block.mat'));
        load(fullfile(expDir.folder,expDir.name));
        %% get photodiode time
        % active task
        win = [0,5000];
        allPD1 = getPhotodiodeTime(session_root,win);
        % ntrials = numel(block.events.feedbackValues);
        ntrials = numel(block.events.endTrialValues);
        allPD1 = allPD1(1:ntrials);
        norepeatValues = block.events.repeatNumValues(1:ntrials);
        norepeat_indx = (norepeatValues==1);
        allPD2 = allPD1(norepeat_indx);
        ntrials_all = ntrials_all+numel(allPD2);
        %% load atlas registration
        load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
        sizeTemplate = [1320,1140];
        Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
        %
        for ii = 1:numel(allPD2)
            stimOn(ii,1) = find(t-allPD2(ii)>0, 1, 'first');
        end
        frames1 = [-70:70];
        frames  = repmat(frames1,[numel(stimOn),1]);
        frames_stimOn = frames+stimOn;
        %% load spirals
        load(fullfile(spiral_folder,'spirals_grouping',[fname '_spirals_group_fftn.mat']));    
        spiral_duration = cellfun(@(x) size(x,1), archiveCell);
        indx2 = (spiral_duration>=2);
        groupedCells = archiveCell(indx2);
        pwAll = cell2mat(groupedCells); 
        %% organize sprials by frame
        [pwAll(:,1),pwAll(:,2)] = transformPointsForward(...
            tform,pwAll(:,1),pwAll(:,2)); 
        pwAll(:,1:2) = round(pwAll(:,1:2));
        [lia,locb] = ismember(pwAll(:,1:2),brain_index,'rows');
        pwAll = pwAll(lia,:);
        spiral_cell = cell(size(frames_stimOn));
        for j = 1:size(frames_stimOn,1)
            for k = 1:size(frames_stimOn,2)
                clear indx1 current_spirals
                current_frame = frames_stimOn(j,k);
                indx1 = find(pwAll(:,5) == current_frame);
                current_spirals = pwAll(indx1,:);
                spiral_cell{j,k} = current_spirals;
            end
        end    
        %%
        spiral_all = cat(1,spiral_all,spiral_cell);
    end
    %%
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
    %%
    correct_index_stimL = (T_all.label == "correct" & T_all.left_contrast >0);
    correct_index_stimR = (T_all.label == "correct" & T_all.right_contrast >0);
    incorrect_index_stimL = (T_all.label == "incorrect" & T_all.left_contrast >0);
    incorrect_index_stimR = (T_all.label == "incorrect" & T_all.right_contrast >0);
    miss_index = (T_all.label == "miss");

    spiral_correct_stimL = spiral_all(correct_index_stimL ,:);
    spiral_correct_stimR = spiral_all(correct_index_stimR ,:);
    spiral_incorrect_stimL = spiral_all(incorrect_index_stimL ,:);
    spiral_incorrect_stimR = spiral_all(incorrect_index_stimR ,:);
    spiral_miss = spiral_all(miss_index,:);

    [spiral_correct_L] = getConcatTrials(spiral_correct_stimL);
    [spiral_correct_R] = getConcatTrials(spiral_correct_stimR);
    [spiral_incorrect_L] = getConcatTrials(spiral_incorrect_stimL);
    [spiral_incorrect_R] = getConcatTrials(spiral_incorrect_stimR);
    [spiral_miss2] = getConcatTrials(spiral_miss);

    correct_L_trialN = sum(correct_index_stimL);
    correct_R_trialN = sum(correct_index_stimR);
    incorrect_L_trialN = sum(incorrect_index_stimL);
    incorrect_R_trialN = sum(incorrect_index_stimR);
    miss_trialN = sum(miss_index);
    %%
    save(fullfile(save_folder,[mn '_spirals_task_sort']),...
        'T_all','spiral_all',...
        'spiral_correct_L','spiral_correct_R',...
        'spiral_incorrect_L','spiral_incorrect_R','spiral_miss2',...
        'correct_L_trialN','correct_R_trialN',...
        'incorrect_L_trialN','incorrect_R_trialN','miss_trialN');
end