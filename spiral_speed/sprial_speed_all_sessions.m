% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
BW = logical(projectedAtlas1);
%% mask for whole brain 
[row,col] = find(BW);
brain_index = [col,row];
%%
lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
Fs = 35; % frame sampling rate
scale1 = 8;
padding = 30;
halfpadding = padding/2;
grid_x = [250,350];
grid_y = [350,550];
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions3.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
for kk = 1:15
    clear filteredSpirals2 filteredSpirals3 tracePhase_raw tracePhase_padded...
        angle_offset_all distance_offset_all  
    %%
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    [U,V,t,mimg] = get_wf_svd(serverRoot);
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %%
    folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
    fname = [mn '_' tdb '_' num2str(en)];
    dfolder = fullfile(folder,mn,td,num2str(en));
    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));
    spiralGroupName = dir(fullfile(dfolder,'*spirals_group_fftn.mat')).name;
    load(fullfile(dfolder,spiralGroupName));
    %%
    clear indx2 groupedCells
    indx2 = cellfun(@(x) size(x,1), archiveCell);
    groupedCells = archiveCell(indx2>=2);
    filteredSpirals2 = cell2mat(groupedCells);
    %%
    index1 = (filteredSpirals2(:,2)>grid_x(1) & filteredSpirals2(:,2)<grid_x(2)...
        &filteredSpirals2(:,1)>grid_y(1) & filteredSpirals2(:,1)<grid_y(2));
    filteredSpirals3 = filteredSpirals2(index1,:);
    %%
    [~,~,tracePhase_raw] = spiralPhaseMap(U(1:scale1:end,1:scale1:end,1:50),dV(1:50,:),t,lowpass);
    tracePhase_raw = permute(tracePhase_raw,[2,3,1]);
    [tracePhase_padded] = padZeros(tracePhase_raw,halfpadding);
    %%
    [angle_offset_all,distance_offset_all] = get_anglular_speed(filteredSpirals3,tracePhase_padded,scale1,halfpadding);
    distance_offset_all1 = cellfun(@(x) mean(x(:,end),'omitnan'),distance_offset_all,'UniformOutput',false);
    distance_offset_all2 = vertcat(distance_offset_all1{:});
    angle_offset_all1 = cellfun(@(x) mean(x(:,end),'omitnan'),angle_offset_all,'UniformOutput',false);
    angle_offset_all2 = vertcat(angle_offset_all1{:});
    mean_angle_offset = cellfun(@(x) mean(x,1),angle_offset_all,'UniformOutput',false);
    mean_angle_all = [mean_angle_offset{:}];
    mean_distance_offset = cellfun(@(x) mean(x,1),distance_offset_all,'UniformOutput',false);
    save(fname,'angle_offset_all','distance_offset_all');
    %%
%     h1 = plot_speed(mean_angle_offset,mean_distance_offset,angle_offset_all2,distance_offset_all2);
%     print(h1,fname, '-dpdf', '-bestfit', '-painters');
%     close all;
end