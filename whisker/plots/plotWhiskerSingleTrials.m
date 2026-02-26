function hs15a = plotWhiskerSingleTrials(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
BW2 = BW(1:2:end,1:2:end);
%%
load(fullfile(data_folder,'whisker','single_trials','single_trial_maps2.mat'));
wf2 = imresize(wf1,[660,570]);
%%
% [hs15a,meanTrace2,tracePhase] = plotMeanTracePhase6(wf2,BW2,maskPath,st,atlas1);
[hs15a] = plotExampleTracePhase6(wf2,BW2,maskPath,st,atlas1);
print(hs15a,fullfile(save_folder,['figS15a_ZYE_0092_example_trials2.pdf']),'-dpdf', '-bestfit', '-painters');