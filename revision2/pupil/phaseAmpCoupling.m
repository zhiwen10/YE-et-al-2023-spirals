function [phase_centers,mean_amp_per_bin_all2] = phaseAmpCoupling(T,data_folder,pixel,behavior, n_bins,low_freq_band,high_freq_band)
%% Parameters
fs = 35;                    % Sampling frequency
mean_amp_per_bin_all = [];
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder
    fname = [mn '_' tdb '_' num2str(en)];
    subfolder = [mn '_' tdb '_' num2str(en)];
    serverRoot = expPath(mn, td, en);
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    % registration
    load(fullfile(data_folder,'spirals','rf_tform',...
        [fname '_tform.mat']));  
    Utransformed = imwarp(U,tform,'OutputView',imref2d([1320,1140]));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d([1320,1140]));
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
    %% Design filters for frequency bands
    % Low frequency (phase) - 4th order Butterworth
    [b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
    signal_low = filtfilt(b_low, a_low, signal2);    
    phase_low = angle(hilbert(signal_low));
    [b_high, a_high] = butter(4, high_freq_band/(fs/2), 'bandpass');
    phase_bins = linspace(-pi, pi, n_bins+1);
    for j = 1:7
        trace = squeeze(Utransformed(pixel(j,1),pixel(j,2),1:50))'*V(1:50,:);
        signal = double(trace)'./mimgtransformed(pixel(j,1),pixel(j,2));
        signal_high = filtfilt(b_high, a_high, signal);
        amplitude_high = abs(hilbert(signal_high));
        mean_amp_per_bin = zeros(n_bins, 1);
        for i = 1:n_bins
            phase_mask = (phase_low >= phase_bins(i)) & (phase_low < phase_bins(i+1));
            mean_amp_per_bin(i) = mean(amplitude_high(phase_mask));
        end
        mean_amp_per_bin_all(:,kk,j) = mean_amp_per_bin;
    end
end
mean_amp_per_bin = mean(mean_amp_per_bin_all,1);
mean_amp_per_bin_all2 = mean_amp_per_bin_all-mean_amp_per_bin; % subtract baseline
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;