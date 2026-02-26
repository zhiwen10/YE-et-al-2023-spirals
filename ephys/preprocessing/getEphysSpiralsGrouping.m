%%
function getEphysSpiralsGrouping(T,data_folder,save_folder)
for kk = 1:size(T,1)
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' opsnum2str(en)];
    load(fullfile(data_folder,'ephys','spirals_raw',[fname '_spirals.mat']));
    filteredSpirals = pwAll(pwAll(:,3)>=40,:);                             % only use sprials with radius >40 pixels, based on 3d-fft
    %% temporal grouping
    filteredSpirals =unique(filteredSpirals, 'rows');                      % get rid of duplication, in case any
    filteredSpirals = sortrows(filteredSpirals,5);                         % sort based on frame number
    [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);       % main grouping algorithm
    save(fullfile(save_folder,[fname '_spirals_group_fftn.mat']),...
        'archiveCell');                                                    % save all archived Cells 
end

