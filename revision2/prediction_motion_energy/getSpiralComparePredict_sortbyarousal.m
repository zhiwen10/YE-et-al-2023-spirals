function getSpiralComparePredict_sortbyarousal(T,data_folder,save_folder)
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
%%
raw_folder = fullfile(data_folder,'ephys','spirals_raw_fftn');
predict_folder = fullfile(data_folder,'ephys','spirals_predict');
reg_folder = fullfile(data_folder,'ephys','rf_tform');
amp_folder = fullfile(data_folder,'ephys','amplitude');
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
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
T = T(logical(T.facevideo),:);
T([14,19],:) = [];
%%
low_freq_band = [0.05, 0.5];   
Fs = 35;
spiral_left_match_all = {};
spiral_right_match_all = {};
ratio_all = [];
ratio =[];
for kk = 1:size(T,1)
    clear spiralsT spirals_filt spirals_prediction_left ...
        spirals_prediction_right spiralsT1 sprials_left spirals_right
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    t = readNPY(fullfile(ops.session_root,...
        'svdTemporalComponents_corr.timestamps.npy'));
    load(fullfile(raw_folder,[fname '_spirals_group_fftn.mat']));
    load(fullfile(reg_folder,[fname '_tform.mat']));
    load(fullfile(amp_folder,[fname '_amp.mat']));
    spiral_length = cellfun(@(x) size(x,1), archiveCell);
    spiral_sequence = archiveCell(spiral_length>=2);
    spirals_filt= cell2mat(spiral_sequence);
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,spirals_filt(:,1),spirals_filt(:,2));
    spirals_filt(:,1:2) = round(spiralsT); 
    %%
    load(fullfile(data_folder,'revision2','prediction_motion_energy',...
    [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    if numel(image_energy2)<size(t,1)
        image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
    elseif numel(image_energy2)>size(t,1)
        image_energy2 = image_energy2(1:numel(t));
    end 
    signal2 = image_energy2;
    [f1,f2] = butter(4, low_freq_band/(Fs/2), 'bandpass');
    meanTrace_low = filtfilt(f1,f2,signal2);
    traceHilbert_low =hilbert(meanTrace_low);
    tracePhase_low = angle(traceHilbert_low);
    %%
    clear pwAll
    load(fullfile(predict_folder,[fname '_spirals_predicted.mat']));
    states = {'on','off'};
    for i = 1:2
        clear spirals_filt1 a b 
        state = states{i};
        if strcmp(state,'on')
            index = find(tracePhase_low >= -pi/2 & tracePhase_low < pi/2);
        elseif strcmp(state,'off')
            index = find(tracePhase_low < -pi/2 | tracePhase_low >= pi/2);
        end
        [a,b] = ismember(spirals_filt(:,5),index);
        spirals_filt1 = spirals_filt(a,:);
        [spiral_left_match,spiral_right_match] =  matchSpirals(spirals_filt1,...
            brain_index_left,brain_index_right,pwAll,tform,traceAmp);
        spiral_left_match_all{kk,i} = spiral_left_match;                        
        spiral_right_match_all{kk,i} = spiral_right_match;       
    end
end
save(fullfile(save_folder,'spiral_compare_sessions_arousal.mat'),...
    'spiral_left_match_all','spiral_right_match_all');
