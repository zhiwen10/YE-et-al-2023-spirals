function spirals_sort = phaseSpiralHistogram(T,data_folder,behavior,n_bins,brain_index,low_freq_band)
%%
Fs = 35;
phase_bins = linspace(-pi, pi, n_bins+1);
spirals_sort =  cell(size(T,1),n_bins);
%%
for kk = 1:size(T,1)                                                 % plot 6 session in one figure, to avoid crowding
    %% session info
    clear spiralsT filteredSpirals trace
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    serverRoot = expPath(mn, td, en);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    %%
    if strcmp(behavior,'pupil')
        tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
        pupil_mean = readNPY(tFile1); 
        pupil_mean2 = pupil_mean(1:2:end);
        if numel(pupil_mean2)<size(t,1)
            pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
        elseif numel(pupil_mean2)>size(t,1)
            pupil_mean2 = pupil_mean2(1:numel(t));
        end 
        signal2 = pupil_mean2;
    elseif strcmp(behavior,'face')
        load(fullfile(data_folder,'spirals','spirals_index',...
        [fname '_motion_energy.mat']));
        image_energy2(isnan(image_energy2)) = 0;
        if numel(image_energy2)<size(t,1)
            image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
        elseif numel(image_energy2)>size(t,1)
            image_energy2 = image_energy2(1:numel(t));
        end 
        signal2 = image_energy2;
    end
%%
    [f1,f2] = butter(4, low_freq_band/(Fs/2), 'bandpass');
    meanTrace_low = filtfilt(f1,f2,signal2);
    traceHilbert_low =hilbert(meanTrace_low);
    tracePhase_low = angle(traceHilbert_low);
    traceAmp_low = abs(traceHilbert_low);
    %% only use spiral sequences with at least 2 consecutive frames
    load(fullfile(data_folder,'spirals','spirals_grouping',...
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