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
low_freq_band = [0.1, 0.5]; % Phase-providing frequency (slow oscillation)
high_freq_band = [2, 8];    % Amplitude-modulated frequency (fast oscillation)
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
    trace = squeeze(Utransformed(pixel(1,1),pixel(1,2),1:50))'*V(1:50,:);
    %%
    % [tAll, tUp, tDown] = getTLtimesDigital(mn, td, en, 'cameraTrigger');
    tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
    pupil_mean = readNPY(tFile1); 
    % pupilsize = [tUp,pupil_mean];
    %%
    pupil_mean2 = pupil_mean(1:2:end);
    if numel(pupil_mean2)<size(t,1)
        pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
    elseif numel(pupil_mean2)>size(t,1)
        pupil_mean2 = pupil_mean2(1:numel(t));
    end 
    %%
    load(fullfile(data_folder,'spirals','spirals_index',...
        [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    %% Assume data is loaded as: signal = your_timeseries_data;
    signal = double(trace)';
    %%
%     f = 0.4;
%     signal = sin(2*pi*f*t);
%     signal2 = wgn(numel(t),1,0);
%     figure;
%     plot(t,signal)
    %% Method 1: Modulation Index (Tort et al., 2010) - Most common approach

    % Design filters for frequency bands
    % Low frequency (phase) - 4th order Butterworth
    [b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
    signal_low = filtfilt(b_low, a_low, signal);

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
    MI = (H_max - H) / H_max; % Modulation Index (0 = no coupling, 1 = perfect coupling)

    %% Method 2: Phase-Locking Value (PLV) between phase and amplitude
    % Convert amplitude to phase (instantaneous phase of amplitude envelope)
    phase_amplitude = angle(hilbert(amplitude_high));
    PLV = abs(mean(exp(1i * (phase_low - phase_amplitude))));

    %% Method 3: Circular-Linear Correlation
    % Correlate circular phase with linear amplitude
    [r_circ, p_circ] = circ_corrcc(phase_low, amplitude_high);

    %% Method 4: Mean Vector Length (MVL)
    % Alternative measure of phase-amplitude coupling
    complex_valued = amplitude_high .* exp(1i * phase_low);
    MVL = abs(mean(complex_valued)) / mean(amplitude_high);

    %% Statistical Testing using Surrogate Data
    n_surrogates = 1000;
    MI_surrogates = zeros(n_surrogates, 1);
    PLV_surrogates = zeros(n_surrogates, 1);
    MVL_surrogates = zeros(n_surrogates, 1);

    fprintf('Computing surrogate statistics for PAC significance testing...\n');
    for s = 1:n_surrogates
        % Time-shift surrogate (randomly shift one signal)
        shift_amount = randi([fs, length(signal) - fs]); % Shift by 1-N seconds
        amplitude_high_surr = circshift(amplitude_high, shift_amount);

        % Compute MI for surrogate
        mean_amp_surr = zeros(n_bins, 1);
        for i = 1:n_bins
            phase_mask = (phase_low >= phase_bins(i)) & (phase_low < phase_bins(i+1));
            mean_amp_surr(i) = mean(amplitude_high_surr(phase_mask));
        end
        P_surr = mean_amp_surr / sum(mean_amp_surr);
        H_surr = -sum(P_surr .* log(P_surr + eps));
        MI_surrogates(s) = (H_max - H_surr) / H_max;

        % PLV for surrogate
        phase_amplitude_surr = angle(hilbert(amplitude_high_surr));
        PLV_surrogates(s) = abs(mean(exp(1i * (phase_low - phase_amplitude_surr))));

        % MVL for surrogate
        complex_surr = amplitude_high_surr .* exp(1i * phase_low);
        MVL_surrogates(s) = abs(mean(complex_surr)) / mean(amplitude_high_surr);
    end

    % Calculate p-values
    p_MI = sum(MI_surrogates >= MI) / n_surrogates;
    p_PLV = sum(PLV_surrogates >= PLV) / n_surrogates;
    p_MVL = sum(MVL_surrogates >= MVL) / n_surrogates;

    %% Method 5: Comodulogram - Comprehensive frequency-by-frequency analysis
    % Test PAC across different frequency combinations
    phase_freqs = 0.1:0.1:1.0;    % Phase frequencies
    amp_freqs = 1.5:0.5:8.0;      % Amplitude frequencies
    comodulogram = zeros(length(phase_freqs), length(amp_freqs));

    fprintf('Computing comodulogram (this may take several minutes)...\n');
    for i = 1:length(phase_freqs)
        for j = 1:length(amp_freqs)
            % Filter for current frequency pair
            [b_p, a_p] = butter(4, [max(0.05, phase_freqs(i)-0.05), phase_freqs(i)+0.05]/(fs/2), 'bandpass');
            [b_a, a_a] = butter(4, [amp_freqs(j)-0.25, min(amp_freqs(j)+0.25, fs/2-0.1)]/(fs/2), 'bandpass');

            try
                sig_phase = filtfilt(b_p, a_p, signal);
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

    %% Visualization
    figure('Position', [100, 100, 1400, 1000]);

    % Plot 1: Phase-amplitude relationship
    subplot(2, 3, 1);
    phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
    bar(phase_centers * 180/pi, mean_amp_per_bin);
    xlabel('Phase of Low Freq (degrees)');
    ylabel('Mean High Freq Amplitude');
    title(sprintf('Phase-Amplitude Coupling\nMI = %.4f (p = %.3f)', MI, p_MI));
    grid on;

    % Plot 2: Polar plot of phase-amplitude coupling
    subplot(2, 3, 2);
    polarplot([phase_centers, phase_centers(1)], [mean_amp_per_bin; mean_amp_per_bin(1)], 'bo-', 'LineWidth', 2);
    title('Polar View of PAC');

    % Plot 3: Time series with overlays
    subplot(2, 3, 3);
    t_display = (1:min(10*fs, length(signal))) / fs; % Show first 10 seconds
    plot(t_display, signal_low(1:length(t_display)), 'b-', 'LineWidth', 1);
    hold on;
    plot(t_display, amplitude_high(1:length(t_display))/max(amplitude_high)*max(signal_low), 'r-', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Low Freq Signal & High Freq Envelope');
    legend('Low Freq Phase', 'High Freq Amplitude', 'Location', 'best');
    grid on;

    % Plot 4: Comodulogram
    subplot(2, 3, 4);
    imagesc(amp_freqs, phase_freqs, comodulogram);
    colorbar;
    xlabel('Amplitude Frequency (Hz)');
    ylabel('Phase Frequency (Hz)');
    title('Comodulogram (Phase-Amp Coupling)');
    axis xy;

    % Plot 5: Statistical comparison
    subplot(2, 3, 5);
    measures = [MI, PLV, MVL];
    p_values = [p_MI, p_PLV, p_MVL];
    surrogate_means = [mean(MI_surrogates), mean(PLV_surrogates), mean(MVL_surrogates)];

    x_pos = 1:3;
    bar(x_pos, measures, 'b');
    hold on;
    bar(x_pos, surrogate_means, 'r');
    set(gca, 'XTickLabel', {'MI', 'PLV', 'MVL'});
    ylabel('Coupling Strength');
    title('PAC Measures vs Surrogate');
    legend('Observed', 'Surrogate Mean', 'Location', 'best');

    % Add p-values as text
    for i = 1:3
        text(i, measures(i) + 0.01, sprintf('p=%.3f', p_values(i)), ...
             'HorizontalAlignment', 'center');
    end

    % Plot 6: Surrogate distributions
    subplot(2, 3, 6);
    histogram(MI_surrogates, 50, 'Normalization', 'pdf', 'FaceAlpha', 0.7);
    hold on;
    xline(MI, 'r-', 'LineWidth', 3, 'Label', sprintf('Observed MI = %.4f', MI));
    xlabel('Modulation Index');
    ylabel('Probability Density');
    title('Surrogate Distribution');
    grid on;

    %% Print Results Summary
    fprintf('\n=== Phase-Amplitude Coupling Results ===\n');
    fprintf('Frequency bands tested:\n');
    fprintf('  Phase band: %.1f - %.1f Hz\n', low_freq_band(1), low_freq_band(2));
    fprintf('  Amplitude band: %.1f - %.1f Hz\n', high_freq_band(1), high_freq_band(2));
    fprintf('\nCoupling measures:\n');
    fprintf('  Modulation Index (MI): %.4f (p = %.3f)\n', MI, p_MI);
    fprintf('  Phase-Locking Value (PLV): %.4f (p = %.3f)\n', PLV, p_PLV);
    fprintf('  Mean Vector Length (MVL): %.4f (p = %.3f)\n', MVL, p_MVL);
    if ~isempty(r_circ)
        fprintf('  Circular-Linear Correlation: %.4f (p = %.3f)\n', r_circ, p_circ);
    end

    fprintf('\nInterpretation:\n');
    if p_MI < 0.05
        fprintf('*** SIGNIFICANT phase-amplitude coupling detected! ***\n');
        fprintf('The phase of %.1f-%.1f Hz oscillations modulates\n', low_freq_band(1), low_freq_band(2));
        fprintf('the amplitude of %.1f-%.1f Hz oscillations.\n', high_freq_band(1), high_freq_band(2));
    else
        fprintf('No significant phase-amplitude coupling detected.\n');
    end

    % Find peak in comodulogram
    [max_pac, max_idx] = max(comodulogram(:));
    [max_i, max_j] = ind2sub(size(comodulogram), max_idx);
    fprintf('\nStrongest PAC in comodulogram:\n');
    fprintf('  Phase freq: %.1f Hz, Amp freq: %.1f Hz, MI = %.4f\n', ...
            phase_freqs(max_i), amp_freqs(max_j), max_pac);
end