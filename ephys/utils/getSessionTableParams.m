%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new2.csv'));
%%
BW = logical(projectedAtlas1);
for kk = 1:35
    % load widefield and spiking data
    ops = get_session_info2(T,kk,data_folder);
    filename = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_mean_explained_var.mat'];
    load(fullfile(data_folder,'ephys','varianceExplainedRegistered',filename));
    var_mean_reg(~BW) = nan;
    max_var(kk,1) = squeeze(max(var_mean_reg,[],[1,2]));
    mean_var(kk,1) = squeeze(mean(var_mean_reg,[1,2],'omitnan'));
    %%
    t = readNPY(fullfile(ops.session_root,'svdTemporalComponents_corr.timestamps.npy'));
    sampleSize(kk,1) = numel(t);
end
sampleT = floor(sampleSize./35);
%%
T1 = T(:,1:10);
T1.sampleSize = sampleSize;
T1.sampleT = sampleT;
T1.meanVar = mean_var;
T1.maxVar = max_var;
%%
writetable(T1,fullfile(data_folder,'tables','spirals_ephys_sessions_new2.csv'));