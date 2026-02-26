function compareSpiralMaps(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data'; 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
raw_folder = fullfile(data_folder,'spirals','spirals_grouping');
predict_folder = fullfile(data_folder,'revision2','spirals_wf_predict');
reg_folder = fullfile(data_folder,'spirals','rf_tform');
%% get index for left and right hemisphere 
scale = 1;
BW = logical(projectedAtlas1);
frameLim = 0;
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 500 400]);
for kk = 1
% for kk = 1:size(T,1)
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
    frame_all = numel(t);
    
    ax1 = subplot(1,2,(kk-1)*2+1);
    [unique_spirals,unique_spirals_unit] = getSpiralDensity3(archiveCell,tform,brain_index,frameLim,frame_all);
    [ax1,cb1]= plotDesnity(ax1,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
    caxis([0,1]);
    %% load spirals from prediction
    clear archiveCell spiral_length spiral_sequence spirals_filt2 spiralsT2
    load(fullfile(data_folder,'revision2','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    spiral_length = cellfun(@(x) size(x,1), archiveCell);
    spiral_sequence = archiveCell(spiral_length>=2);
    spirals_filt2 = cell2mat(spiral_sequence);
    spirals_filt2(spirals_filt2(:,4)==0,:) = -1;
    [spiralsT2(:,1),spiralsT2(:,2)] = transformPointsForward(...
        tform,spirals_filt2(:,1),spirals_filt2(:,2));
    spirals_filt2(:,1:2) = round(spiralsT2);   
    
    ax2 = subplot(1,2,(kk-1)*2+2);
    [unique_spirals,unique_spirals_unit] = getSpiralDensity3(archiveCell,tform,brain_index,frameLim,frame_all);
    [ax1,cb1]= plotDesnity(ax2,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
    caxis([0,1]);
end

