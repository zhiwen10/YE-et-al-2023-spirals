function getSpiralComparePredictWF(T,data_folder,save_folder)
% for each sprials in the raw spiral list, find a (and only one) spiral 
% in the nearby 100 pixel distance in the predicted frame. 
% If exist, then save both raw and preidcted spiral as a row. 

% In the final saved cell structure (left and right): 
% each cell is a single session
% Within each cell, each row is a spiral that has a nearby spiral in
% prediction

% Col  1       2        3      4          5      
%    raw_x    raw_y  radius  direction  frameN

% Col  6       7        8      9          10      
%    pdt_x    pdt_y  radius  direction  frameN

% Col  11           12
%    match   frame 2-8Hz amp
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
raw_folder = fullfile(data_folder,'spirals','spirals_grouping');
predict_folder = fullfile(data_folder,'revision2','spirals_wf_predict');
reg_folder = fullfile(data_folder,'spirals','rf_tform');
%% get index for left and right hemisphere 
scale = 1;
BW = logical(projectedAtlas1);
BW_right = BW; BW_right(:,1:size(projectedAtlas1,2)/2) = 0; 
BW_left = BW; BW_left(:,size(projectedAtlas1,2)/2+1:end) = 0; 
[rowL,colL] = find(BW_left);
brain_index_left = [colL,rowL];
[rowR,colR] = find(BW_right);
brain_index_right = [colR,rowR];
%%
spiral_left_match_all = {};
spiral_right_match_all = {};
ratio_all = [];
ratio =[]; 
for kk = 1:size(T,1)
    clear spiralsT spirals_filt spirals_prediction_left ...
        spirals_prediction_right spiralsT1 sprials_left spirals_right
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder   
    subfolder = [mn '_' tdb '_' num2str(en)];
    fname = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    t = readNPY(fullfile(session_root,...
        'svdTemporalComponents_corr.timestamps.npy'));
    load(fullfile(raw_folder,[fname '_spirals_group_fftn.mat']));
    load(fullfile(reg_folder,[fname '_tform.mat']));
    spiral_length = cellfun(@(x) size(x,1), archiveCell);
    spiral_sequence = archiveCell(spiral_length>=2);
    spirals_filt= cell2mat(spiral_sequence);
    spirals_filt(spirals_filt(:,4)==0,:) = -1;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,spirals_filt(:,1),spirals_filt(:,2));
    spirals_filt(:,1:2) = round(spiralsT); 
    %% search for match in left and right hemispehre seperately
    [liaL,locbL] = ismember(spirals_filt(:,1:2),brain_index_left,'rows');
    [liaR,locbR] = ismember(spirals_filt(:,1:2),brain_index_right,'rows');
    spirals_left = spirals_filt(liaL,:);
    spirals_right = spirals_filt(liaR,:);
    %% load spirals from prediction
    clear pwAll
    T2 = readtable(fullfile(data_folder,'revision2','spirals_wf_predict',...
        [fname '_spirals_all.csv']));
    pwAll = table2array(T2);
    [spiralsT1(:,1),spiralsT1(:,2)] = transformPointsForward(...
        tform,pwAll(:,1),pwAll(:,2));
    pwAll(:,1:2) = round(spiralsT1); 
    [liaL,locbL] = ismember(pwAll(:,1:2),brain_index_left,'rows');
    spirals_prediction_left = pwAll(liaL,:);
    [liaR,locbR] = ismember(pwAll(:,1:2),brain_index_right,'rows');
    spirals_prediction_right = pwAll(liaR,:);
    %% find spirals within 100 pixels distance from each raw spiral
    % only care about the situation when only one spiral is within this
    % distance
    clear match spiral_p1 spiral_p2 index1 sprial_index
    spiral_left_match = [];                                                % check left hemisphere first 
    for i = 1:size(spirals_left)
        frame = spirals_left(i,5);
        spiralx = spirals_left(i,1); spiraly = spirals_left(i,2);
        
        spiral_index = find(spirals_prediction_left(:,5)==frame);
        spiral_p1 = spirals_prediction_left(spiral_index,:);
        index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
            spiral_p1(:,1)<= spiralx+100 & ...
            spiral_p1(:,2)>= spiraly-100 & ...
            spiral_p1(:,2)<= spiraly+100);
        spiral_p2 = spiral_p1(index1,:);
        if size(spiral_p2,1) == 1
            spiral_temp = [spirals_left(i,1:5), spiral_p2];
            spiral_left_match = [spiral_left_match;spiral_temp];
        end
    end
    spiral_left_match(:,11) = 0;
    match = spiral_left_match(:,9)-spiral_left_match(:,4);                 % check if spiral directions are matching in raw and predicted
    match = not(match);
    spiral_left_match(:,11) = match;
    spiral_left_match_all{kk} = spiral_left_match;                         % save all assessed sprials in a session to the cell structure 
    %%
    clear match spiral_p1 spiral_p2 index1 sprial_index
    spiral_right_match = [];                                               % now check right hemisphere 
    for i = 1:size(spirals_right)
        frame = spirals_right(i,5);
        spiralx = spirals_right(i,1); spiraly = spirals_right(i,2);
        
        spiral_index = find(spirals_prediction_right(:,5)==frame);
        spiral_p1 = spirals_prediction_right(spiral_index,:);
        index1 = find(spiral_p1(:,1)>= spiralx-100 & ...
            spiral_p1(:,1)<= spiralx+100 & ...
            spiral_p1(:,2)>= spiraly-100 & ...
            spiral_p1(:,2)<= spiraly+100);
        spiral_p2 = spiral_p1(index1,:);
        if size(spiral_p2,1) == 1
            spiral_temp = [spirals_right(i,1:5), spiral_p2];
            spiral_right_match = [spiral_right_match;spiral_temp];
        end
    end
    spiral_right_match(:,11) = 0;
    match = spiral_right_match(:,9)-spiral_right_match(:,4);               % check if spiral directions are matching in raw and predicted
    match = not(match);
    spiral_right_match(:,11) = match;
    spiral_right_match_all{kk} = spiral_right_match;                       % save all assessed sprials in a session to the cell structure 
end
%%
save(fullfile(save_folder,'spiral_compare_sessions_neighbor.mat'),...
    'spiral_left_match_all','spiral_right_match_all');
