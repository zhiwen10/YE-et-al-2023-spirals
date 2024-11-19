function getSpiralDurationRatio_freq(T,freq,data_folder,save_folder)
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
for kk = 1:size(T,1)
    clear p pwAll pwAll1 archiveCell
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% get spiral ratio sorted by the spiral sequence length
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_freq','spirals_fftn',...
        freq_folder,[fname '_spirals_group_fftn.mat']));                               % load grouped spiral centers (>40 pixels radius)
    indx2 = cellfun(@(x) size(x,1), archiveCell);
    spirals_total = size(cell2mat(archiveCell),1);
    for duration = 1:50
        indx3 = (indx2==duration);                                         % find all grouped cells that had N sprials in a sequence 
        sprial_seq_temp = archiveCell(indx3);                   
        N(duration,1) = size(cell2mat(sprial_seq_temp),1);                 % count number of all spirals
        N_ratio(duration,1) = N(duration,1)./spirals_total;                % calculate spirals/ all spirals
    end
    %% get spiral ratio sorted by the spiral sequence length in permutation
    N_scramble = [];
    N_ratio_scramble = [];
    for ii = 1:10
        clear indx_scramble
        filename = fullfile(data_folder,'spirals','spirals_freq','spirals_scrambled',...
            freq_folder,fname,[fname '_scramble_group_' num2str(ii)]);
        load(filename);                                                    % load each scrambled spiral group (10 in total)
        indx_scramble = cellfun(@(x) size(x,1), archiveCell);
        for duration = 1:50                                                % find spiral ratio for all 50 length of sequences
            clear indx_scramble_temp
            indx_scramble_temp = (indx_scramble==duration);                % find all grouped cells that had N sprials in a sequence 
            sprial_seq_temp = archiveCell(indx_scramble_temp);
            N_scramble(duration,ii) = size(cell2mat(sprial_seq_temp),1);   % count number of all spirals
            N_ratio_scramble(duration,ii) = N_scramble(duration,ii)./...   % calculate spirals/ all spirals
                spirals_total;
        end
    end
    %%
    save_folder1 = fullfile(save_folder,freq_folder);
    mkdir(save_folder1);
    save(fullfile(save_folder1,fname),...
        'N', 'N_ratio','N_scramble','N_ratio_scramble','spirals_total');
end