%%
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
%%
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
T = T(T.face &T.eye_DLC,:);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% Parameters
fs = 35;                    % Sampling frequency
% Frequency bands
low_freq_band = [0.1, 0.3]; % Phase-providing frequency (slow oscillation)
high_freq_band = [2, 8];    % Amplitude-modulated frequency (fast oscillation)
% high_freq_band = [0.5, 2];    % Amplitude-modulated frequency (fast oscillation)

phase_freqs = 0.1:0.1:1;    % Phase frequencies
amp_freqs = 2:1:8;      % Amplitude frequencies
comodulogram = [];
comodulogram_all = [];
%%
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
    %% registration
    load(fullfile(data_folder,'spirals','rf_tform',...
        [fname '_tform.mat']));  
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    %%
    trace = squeeze(Utransformed(pixel(3,1),pixel(3,2),1:50))'*V(1:50,:);
    %%
    tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
    load(fullfile(data_folder,'spirals','spirals_index',...
    [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    pupil_mean = readNPY(tFile1); 
    pupil_mean2 = pupil_mean(1:2:end);
    if numel(pupil_mean2)<size(t,1)
        pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
    elseif numel(pupil_mean2)>size(t,1)
        pupil_mean2 = pupil_mean2(1:numel(t));
    end 
    if numel(image_energy2)<size(t,1)
        image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
    elseif numel(image_energy2)>size(t,1)
        image_energy2 = image_energy2(1:numel(t));
    end 
    %% Assume data is loaded as: signal = your_timeseries_data;
    signal = double(trace)';
    signal2 = image_energy2;
    %% Method 1: Modulation Index (Tort et al., 2010) - Most common approach
    % Design filters for frequency bands
    % Low frequency (phase) - 4th order Butterworth
    [b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
    signal_low = filtfilt(b_low, a_low, signal2);
    % High frequency (amplitude) - 4th order Butterworth  
    [b_high, a_high] = butter(4, high_freq_band/(fs/2), 'bandpass');
    signal_high = filtfilt(b_high, a_high, signal);
    % Extract phase from low frequency signal using Hilbert transform
    phase_low = angle(hilbert(signal_low));
    % Extract amplitude envelope from high frequency signal
    amplitude_high = abs(hilbert(signal_high));
    % Compute Modulation Index (MI)
    n_bins = 18; % Number of phase bins (20-degree bins)
    phase_bins = linspace(-pi, pi, n_bins+1);
    mean_amp_per_bin = zeros(n_bins, 1);
    for i = 1:n_bins
        phase_mask = (phase_low >= phase_bins(i)) & (phase_low < phase_bins(i+1));
        mean_amp_per_bin(i) = mean(amplitude_high(phase_mask));
    end
    % Calculate Modulation Index
    P = mean_amp_per_bin / sum(mean_amp_per_bin); % Normalize to probability distribution
    H = -sum(P .* log(P + eps)); % Entropy (add eps to avoid log(0))
    H_max = log(n_bins); % Maximum possible entropy
    MI(:,kk) = (H_max - H) / H_max; % Modulation Index (0 = no coupling, 1 = perfect coupling)
    %% Method 5: Comodulogram - Comprehensive frequency-by-frequency analysis
    % Test PAC across different frequency combinations
    % amp_freqs = 0.5:0.2:2;      % Amplitude frequencies
    comodulogram = zeros(length(phase_freqs), length(amp_freqs));

    fprintf('Computing comodulogram (this may take several minutes)...\n');
    for i = 1:length(phase_freqs)
        for j = 1:length(amp_freqs)
            % Filter for current frequency pair
            [b_p, a_p] = butter(4, [max(0.05, phase_freqs(i)-0.05), phase_freqs(i)+0.05]/(fs/2), 'bandpass');
            [b_a, a_a] = butter(4, [amp_freqs(j)-0.5, min(amp_freqs(j)+0.5, fs/2)]/(fs/2), 'bandpass');
            % [b_a, a_a] = butter(4, [amp_freqs(j)-0.1, min(amp_freqs(j)+0.1, fs/2)]/(fs/2), 'bandpass');

            try
                sig_phase = filtfilt(b_p, a_p, signal2);
                sig_amp = filtfilt(b_a, a_a, signal);

                % Compute PAC
                phase_curr = angle(hilbert(sig_phase));
                amp_curr = abs(hilbert(sig_amp));

                % Modulation Index for current frequency pair
                mean_amp_curr = zeros(n_bins, 1);
                for k = 1:n_bins
                    mask = (phase_curr >= phase_bins(k)) & (phase_curr < phase_bins(k+1));
                    if sum(mask) > 0
                        mean_amp_curr(k) = mean(amp_curr(mask));
                    end
                end
                P_curr = mean_amp_curr / sum(mean_amp_curr);
                H_curr = -sum(P_curr .* log(P_curr + eps));
                comodulogram(i, j) = (H_max - H_curr) / H_max;
            catch
                comodulogram(i, j) = 0; % In case of filtering issues
            end
        end
    end
    %%
    comodulogram_all(:,:,kk) = comodulogram;
    % phase_centers_all(:,kk) = phase_centers;
    mean_amp_per_bin_all(:,kk) = mean_amp_per_bin;
end
mean_amp_per_bin_all_sum = sum(mean_amp_per_bin_all,1);
mean_amp_per_bin_all2 = mean_amp_per_bin_all./mean_amp_per_bin_all_sum;
comodulogram_all_mean = mean(comodulogram_all,3);
%%
figure;
% Plot 4: Comodulogram
subplot(1,2,1);
imagesc(amp_freqs, phase_freqs, comodulogram_all_mean);
colorbar;
xlabel('Amplitude Frequency (Hz)');
ylabel('Phase Frequency (Hz)');
title('Comodulogram (Phase-Amp Coupling)');
axis xy;
subplot(1,2,2);
amp_per_bin_mean = mean(mean_amp_per_bin_all2,2);
amp_per_bin_sem = std(mean_amp_per_bin_all2,[],2)./sqrt(size(mean_amp_per_bin_all2,2));
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
errorbar(phase_centers * 180/pi, amp_per_bin_mean,amp_per_bin_sem);
xlabel('Phase of Low Freq (degrees)');
ylabel('Mean High Freq Amplitude');
grid on;
%%
a1 = squeeze(comodulogram_all(2,:,:));
figure;
for i = 1:size(comodulogram_all,3)
    plot(amp_freqs,a1(:,i),'k');
    hold on;
end
plot(amp_freqs,mean(a1,2),'r');