%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new.csv'));
T1 = T(30:35,:);
%%
for kk = 1:6
    %%
    ops = get_session_info2(T1,kk,data_folder);
    reg_name = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x.mat'];
    reg_fullname =fullfile(data_folder,'ephys','rf_tform_4x',reg_name);
    var_name = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_explained_var.mat'];
    var_fullname =fullfile(data_folder,'ephys','varianceExplainedAll',var_name);
    %%
    load(var_fullname);
    load(reg_fullname);
    %%
    var_mean = squeeze(mean(explained_var_all,3));
    sizeTemplate = size(projectedTemplate1);
    var_mean_reg = imwarp(var_mean,tform,'OutputView',imref2d(sizeTemplate));
    filename = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_mean_explained_var'];
    save(fullfile(data_folder,'ephys','varianceExplainedRegistered',filename),'var_mean_reg');
end