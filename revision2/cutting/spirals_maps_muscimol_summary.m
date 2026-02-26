githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'YE-et-al-2023-spirals')));            % paper repository
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
spiral_folder = 'E:\task2\spirals\spirals_grouping';
data_folder = 'E:\task2';
T1 = readtable('cutting_session_list.xlsx');
%%
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
[sessions,ia,ic] = unique(T2.session_id);
T3 = T2(ia,:);
%%
label_all = {'spontaneous','control','cutting'};
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 800]);
for m = 1:numel(sessions)
    %%
    mn = T3.MouseID{m};
    session = sessions(m);    
    %%
    for k = 1:3
        %%
        clear archiveCell filteredSpirals spiralsT lia
        %%
        label = label_all{k};
        Ti = T1(ismember(T1.MouseID, mn) & T1.session_id == session & ismember(T1.label, label),:);

        tda = Ti.date;
        en = Ti.folder;
        tdb = datestr(tda,'yyyymmdd');
        fname = [mn '_' tdb '_' num2str(en)];
        load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));
        load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
        session_root = fullfile(data_folder,'task_svd',fname);
        load(fullfile(session_root,'svdTemporalComponents_corr_timestamps.mat')); 
        frame_all = numel(t);
        frameLim = 0;
        ax1 = subplot(3,3,(m-1)*3+k);
        [unique_spirals,unique_spirals_unit] = getSpiralDensity3(archiveCell,tform,brain_index,frameLim,frame_all);
        [ax1,cb1]= plotDesnity(ax1,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
        caxis([0,2]);
    end
end
%%
print(h1, 'SSp_cutting_summary', '-dpdf', '-bestfit', '-painters');

