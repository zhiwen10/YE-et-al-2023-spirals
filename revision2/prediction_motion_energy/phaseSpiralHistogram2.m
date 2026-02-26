function spirals_sort = phaseSpiralHistogram2(T,data_folder,n_bins,brain_index,low_freq_band)
%%
Fs = 35;
phase_bins = linspace(-pi, pi, n_bins+1);
spirals_sort =  cell(size(T,1),n_bins);
%%
for kk = 1:size(T,1)                                                 % plot 6 session in one figure, to avoid crowding
    %% session info
    clear spiralsT filteredSpirals trace
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    t = readNPY(fullfile(ops.session_root,...
        'svdTemporalComponents_corr.timestamps.npy'));
    %% registration
    load(fullfile(data_folder,'ephys','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    %%
    load(fullfile(data_folder,'revision2','prediction_motion_energy',...
    [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    if numel(image_energy2)<size(t,1)
        image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
    elseif numel(image_energy2)>size(t,1)
        image_energy2 = image_energy2(1:numel(t));
    end 
    signal2 = image_energy2;
%%
    [f1,f2] = butter(4, low_freq_band/(Fs/2), 'bandpass');
    meanTrace_low = filtfilt(f1,f2,signal2);
    traceHilbert_low =hilbert(meanTrace_low);
    tracePhase_low = angle(traceHilbert_low);
    traceAmp_low = abs(traceHilbert_low);
    %% only use spiral sequences with at least 2 consecutive frames
    load(fullfile(data_folder,'ephys','spirals_raw_fftn',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);                                          % sprial length > = 2 frames
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));                  % transform sprials to atlas space
    filteredSpirals(:,1:2) = round(spiralsT);    
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % find spirals within the brain boundry
    filteredSpirals = filteredSpirals(lia,:);   
    for mm = 1:n_bins
        clear index a b       
        index = find(tracePhase_low >= phase_bins(mm) & tracePhase_low < phase_bins(mm+1));
        [a,b] = ismember(filteredSpirals(:,5),index);
        spirals_sort{kk,mm} = filteredSpirals(a,:);
    end  
end