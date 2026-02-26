function getSpiralsGroup_freq(T,data_folder,save_folder)
freq_all = [0.2 0.5;0.5,2;2,8];
for ifreq = 1:3
    freq = freq_all(ifreq,:);
    freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
    save_folder1 = fullfile(save_folder, freq_folder);
    mkdir(save_folder1);
    for kk = 1:size(T,1)
        mn = T.MouseID{kk};
        tda = T.date(kk);
        en = T.folder(kk);
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        %% load raw detected spirals
        fname = [mn '_' tdb '_' num2str(en)];
        load(fullfile(data_folder,'spirals','spirals_freq','raw',...
            freq_folder,[fname '_spirals_all.mat']));
        filteredSpirals = pwAll(pwAll(:,3)>=40,:);                             % only use sprials with radius >40 pixels, based on 3d-fft
        %% temporal grouping
        filteredSpirals =unique(filteredSpirals, 'rows');                      % get rid of duplication, in case any
        filteredSpirals = sortrows(filteredSpirals,5);                         % sort based on frame number
        [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);       % main grouping algorithm
        save(fullfile(save_folder1,[fname '_spirals_group_fftn.mat']),...       % save all archived Cells  
            'archiveCell');
    end
end