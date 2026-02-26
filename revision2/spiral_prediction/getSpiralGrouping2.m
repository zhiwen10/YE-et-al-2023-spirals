function getSpiralGrouping2(T,data_folder,save_folder)
%%
% only use spirals with radius > 40 pixels, 
% as distribution of smaller spirals did not exceed 
% random occurance, based on 3d-fft analysis
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    T2 = readtable(fullfile(data_folder,'revision2','MO_predict','spirals_predict_MO_both',...
        [fname '_spirals_all.csv']));
    pwAll = table2array(T2);
    filteredSpirals = pwAll(pwAll(:,3)>=40,:);                             % only use sprials with radius >40 pixels, based on 3d-fft
    %% temporal grouping
    filteredSpirals =unique(filteredSpirals, 'rows');                      % get rid of duplication, in case any
    filteredSpirals = sortrows(filteredSpirals,5);                         % sort based on frame number
    [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);       % main grouping algorithm
    save(fullfile(save_folder,[fname '_spirals_group_fftn.mat']),...       % save all archived Cells  
        'archiveCell');
end
