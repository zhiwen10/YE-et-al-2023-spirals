function h5b = plotWhiskerMeanMap(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
BW2 = BW(1:2:end,1:2:end);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
%% load mean map across 5 animals
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
load(fullfile(data_folder,'whisker','whisker_mean_maps','whisker_spirals_mean_all.mat'));
wf_mean3 = imresize(wf_mean2,[660,570]);
% load detected sprials                                              
load(fullfile(data_folder,'whisker','whisker_mean_maps','whisker_spirals_group_fftn.mat'));
spirals = cell2mat(archiveCell);
%% tracePhase
wf1 = reshape(wf_mean2,size(wf_mean2,1)*size(wf_mean2,2),size(wf_mean2,3));
meanTrace = wf1 -mean(wf1 ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
meanTrace2 = meanTrace';
meanTrace2 = reshape(meanTrace2,size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));

traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = tracePhase';
tracePhase = reshape(tracePhase, size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));
%%
[h5b,meanTrace2,tracePhase] = plotMeanTracePhase7(wf_mean3,BW2,maskPath,st,atlas1,spirals);
%%
print(h5b,fullfile(save_folder, 'Fig5b_WhiskerEvokedMeanMaps.pdf'),...
    '-dpdf', '-bestfit', '-painters');
