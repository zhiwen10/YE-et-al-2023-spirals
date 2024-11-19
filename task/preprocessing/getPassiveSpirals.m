function getPassiveSpirals(data_folder,save_folder)
%% load atlas brain horizontal projection and outline 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
scale = 1;
[row,col] = find(BW);
brain_index = [col,row];
%%
spiral_folder = fullfile(data_folder,'task','spirals_all');
%%
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
for m = 1:4
    %%
    mn = fnames{m};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "passive",:);
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
        win = [0,5000];
        allPD1 = getPhotodiodeTime2(session_root,block,win);
        ntrials = numel(block.events.endTrialValues);
        allPD2 = allPD1(1:ntrials);
        ntrials_all = ntrials_all+numel(allPD2);
        %% load atlas registration
        load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
        sizeTemplate = [132,1140];
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
        %
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
    win = [0,5000];
    [contrast_all] = getPassiveContrasts(T1,data_folder,win);

    high_index_stimR = (contrast_all(:,2)==1  | contrast_all(:,2)==0.5);
    spiral_high_stimR = spiral_all(high_index_stimR,:);
    low_index_stimR = (contrast_all(:,2)==0.25  | contrast_all(:,2)==0.125);
    spiral_low_stimR = spiral_all(low_index_stimR,:);

    high_index_stimL = (contrast_all(:,1)==1  | contrast_all(:,1)==0.5);
    spiral_high_stimL = spiral_all(high_index_stimL,:);
    low_index_stimL = (contrast_all(:,1)==0.25  | contrast_all(:,1)==0.125);
    spiral_low_stimL = spiral_all(low_index_stimL,:);

    zero_index = (contrast_all(:,1)==0  & contrast_all(:,2)==0);
    spiral_zero = spiral_all(zero_index,:);

    spiral_high_stimR_trialN = sum(high_index_stimR);
    spiral_low_stimR_trialN = sum(low_index_stimR);
    spiral_high_stimL_trialN = sum(high_index_stimL);
    spiral_low_stimL_trialN = sum(low_index_stimL);
    spiral_zero_trialN = sum(zero_index);
    %%
    save(fullfile(save_folder,[mn '_spirals_passive_sort']),'spiral_all',...
        'spiral_high_stimL','spiral_low_stimL',...
        'spiral_high_stimR','spiral_low_stimR','spiral_zero',...
        'spiral_high_stimL_trialN','spiral_low_stimL_trialN',...
        'spiral_high_stimR_trialN','spiral_low_stimR_trialN','spiral_zero_trialN');
end