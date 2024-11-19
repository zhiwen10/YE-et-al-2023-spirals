function getSpiralGroupingScrambled(T,data_folder,save_folder)
%%
% permute the detected sprials based on which frame it is on.
% so that we can keep the structure of spiral frames, 
% but with a permuted spiral frame sequence 
% spirals are then grouped with the same algorithm as raw data,
% and saved as mat files. Each session is permuted 10x
%%
tic;
for kk = 1:size(T,1)
    clear p pwAll pwAll1 archiveCell
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    T2 = readtable(fullfile(data_folder,'spirals','spirals_raw',...
        [fname '_spirals_all.csv']));                                      % load all spirals
    pwAll = table2array(T2);            
    pwAll = pwAll(pwAll(:,3)>=40,:);                                       % only use spirals >=40pixels radius, based on fftn
    for ii = 1:10
        p = randperm(size(pwAll,1));                                       % generate random permutation of index
        frame_number = pwAll(:,5);                                         % column5 is the frameN of the detected sprials
        pwAll1 = pwAll;
        pwAll1(:,5) = frame_number(p);                                     % re-order the frameN for each spiral randomly
        pwAll1 = sortrows(pwAll1,5);
        [archiveCell,test_stats] = getGroupingAlgorithm(pwAll1);           % main grouping algorithm
        filename = fullfile(save_folder,...
            [fname '_scramble_group_' num2str(ii)]);
        save(filename,'archiveCell','p');
        %%
        T1(ii) = toc;
        fprintf('mouseID %g -- iter %g; time elapsed %g s \n',...
            [kk,ii,T1(ii)]);
    end
end
