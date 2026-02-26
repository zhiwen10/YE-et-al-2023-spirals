%%
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
%%
data_folder = 'E:\spiral_data_share\data';   
%%
mn = 'AB_0004';
td = '2021-03-30';
en = 1;
serverRoot = expPath(mn, td, en);
%%
tdb = datestr(td,'yyyymmdd');
subfolder = [mn '_' tdb '_' num2str(en)];
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
%%
figure;
imagesc(mimg)
pixel = [400,420];
trace = squeeze(U(pixel(1),pixel(2),1:50))'*V(1:50,:);
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
%%
figure;
subplot(3,1,1);
plot(t,pupil_mean2);
subplot(3,1,2);
plot(t,image_energy2);
subplot(3,1,3);
plot(t,trace);
%% Frequency Coherence Analysis for Two Time Series
% For 30-minute data sampled at 35 Hz, focusing on 0.1-8 Hz
%% Parameters
fs = 35;                    % Sampling frequency (Hz)
duration = numel(t)/35;         % 30 minutes in seconds
N = fs * duration;          % Total number of samples (63,000)

% Assume your data is loaded as:
% x = pupil_mean2;   % Replace with your actual data
x = image_energy2;
% y = double(trace)';  % Replace with your actual data
y = pupil_mean2;
%% Method 1: Standard mscohere with optimized parameters
% Window parameters for good 0.1-8 Hz resolution
window_length = 2^12;       % 4096 samples = ~117 seconds
                           % Gives freq resolution of 35/4096 = 0.0085 Hz
overlap = window_length/2;  % 50% overlap
nfft = window_length;      % Same as window length

% Hamming window for reduced spectral leakage
window = hamming(window_length);

% Compute coherence
[Cxy, f] = mscohere(x, y, window, overlap, nfft, fs);
%% Method 2: Alternative with different window size for comparison
% Shorter window for better time resolution but coarser frequency resolution
window_length2 = 2^11;      % 2048 samples = ~58 seconds
overlap2 = window_length2/2;
window2 = hamming(window_length2);

[Cxy2, f2] = mscohere(x, y, window2, overlap2, window_length2, fs);

%% Focus on frequency range of interest (0.1-8 Hz)
freq_range = (f >= 0.1) & (f <= 8);
f_roi = f(freq_range);
Cxy_roi = Cxy(freq_range);

freq_range2 = (f2 >= 0.1) & (f2 <= 8);
f_roi2 = f2(freq_range2);
Cxy_roi2 = Cxy2(freq_range2);

%% Statistical significance testing - Multiple approaches

% % Method 1: Theoretical threshold based on degrees of freedom
% num_segments = floor((N - overlap) / (window_length - overlap));
% dof = 2 * num_segments;  % Degrees of freedom
% alpha = 0.05;
% % For coherence, use F-distribution based threshold
% F_crit = finv(1-alpha, 2, dof-2);
% theoretical_threshold = F_crit / (F_crit + dof - 2);
% 
% % Method 2: Surrogate data method (most robust)
% n_surrogates = 1000;
% surrogate_coherences = zeros(n_surrogates, length(f));
% 
% fprintf('Computing surrogate data threshold (this may take a moment)...\n');
% for i = 1:n_surrogates
%     % Phase randomization preserving amplitude spectrum
%     X_fft = fft(x);
%     random_phases = 2*pi*rand(size(X_fft));
%     x_surrogate = real(ifft(abs(X_fft) .* exp(1i*random_phases)));
%     
%     % Compute coherence with surrogate
%     [Cxy_surr, ~] = mscohere(x_surrogate, y, window, overlap, nfft, fs);
%     surrogate_coherences(i, :) = Cxy_surr;
% end
% 
% % 95th percentile as threshold
% surrogate_threshold = prctile(surrogate_coherences, 95, 1);
% 
% % Method 3: Bootstrap confidence intervals
% n_bootstrap = 500;
% bootstrap_coherences = zeros(n_bootstrap, length(f));
% 
% fprintf('Computing bootstrap confidence intervals...\n');
% for i = 1:n_bootstrap
%     % Resample with replacement
%     n_samples = length(x);
%     boot_indices = randi(n_samples, n_samples, 1);
%     x_boot = x(boot_indices);
%     y_boot = y(boot_indices);
%     
%     [Cxy_boot, ~] = mscohere(x_boot, y_boot, window, overlap, nfft, fs);
%     bootstrap_coherences(i, :) = Cxy_boot;
% end
% 
% % Bootstrap 95% confidence intervals
% bootstrap_lower = prctile(bootstrap_coherences, 2.5, 1);
% bootstrap_upper = prctile(bootstrap_coherences, 97.5, 1);
%% Method 4: Simple empirical threshold (often most practical)
% Based on the idea that most frequencies should show low coherence
empirical_threshold = prctile(Cxy, 95); % 95th percentile of all coherence values

fprintf('\n=== Significance Thresholds ===\n');
fprintf('Theoretical (F-test): %.3f\n', theoretical_threshold);
fprintf('Surrogate data (95th percentile): %.3f (mean across frequencies)\n', mean(surrogate_threshold));
fprintf('Bootstrap median threshold: %.3f\n', median(bootstrap_upper));
fprintf('Empirical (95th percentile): %.3f\n', empirical_threshold);

%% Plotting results
figure('Position', [100, 100, 1200, 800]);

% Main coherence plot with multiple thresholds
subplot(2,2,1);
plot(f_roi, Cxy_roi, 'b-', 'LineWidth', 2);
hold on;
yline(empirical_threshold, 'm--', 'LineWidth', 1.5, 'DisplayName', sprintf('Empirical (%.3f)', empirical_threshold));
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Coherence with Multiple Significance Thresholds');
legend('Location', 'best');
grid on;
xlim([0.1, 8]);
ylim([0, 1]);

% Surrogate comparison
subplot(2,2,2);
plot(f_roi2, Cxy_roi2, 'b-', 'LineWidth', 2, 'DisplayName', 'Observed');
hold on;
yline(empirical_threshold, 'm--', 'LineWidth', 1.5, 'DisplayName', sprintf('Empirical (%.3f)', empirical_threshold));
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Observed vs Surrogate Data');
legend('Location', 'best');
grid on;
xlim([0.1, 8]);
ylim([0, 1]);

% Bootstrap confidence intervals
subplot(2,2,3);
semilogx(f_roi, Cxy_roi, 'b-', 'LineWidth', 2, 'DisplayName', 'Observed');
hold on;
yline(empirical_threshold, 'm--', 'LineWidth', 1.5, 'DisplayName', sprintf('Empirical (%.3f)', empirical_threshold));
xlabel('Frequency (Hz) - Log Scale');
ylabel('Coherence');
title('Coherence - Log Frequency Scale');
legend('Location', 'best');
grid on;
xlim([0.1, 10]);
ylim([0, 1]);
xticks([0.1,0.2,0.5,1,2,4,8]);

% Frequency bands analysis
subplot(2,2,4);
% Define physiological frequency bands
bands = struct();
bands.slow_1 = [0.1, 0.5];     % Very slow fluctuations
bands.slow_2 = [0.5, 1.0];     % Slow fluctuations  
bands.low = [1.0, 2.0];        % Low frequency
bands.mid = [2.0, 4.0];        % Mid frequency
bands.high = [4.0, 8.0];       % Higher frequency

band_names = fieldnames(bands);
band_coherence = zeros(length(band_names), 1);

for i = 1:length(band_names)
    band_range = (f >= bands.(band_names{i})(1)) & (f <= bands.(band_names{i})(2));
    band_coherence(i) = mean(Cxy(band_range));
end

bar(band_coherence);
set(gca, 'XTickLabel', {'0.1-0.5', '0.5-1.0', '1.0-2.0', '2.0-4.0', '4.0-8.0'});
xlabel('Frequency Bands (Hz)');
ylabel('Mean Coherence');
title('Average Coherence by Frequency Band');
grid on;
hold on;
yline(confidence_level, 'r--', 'LineWidth', 1.5);

%% Additional analysis: Phase and cross-spectral density
figure('Position', [200, 200, 1000, 600]);

% Compute cross-spectral density and phase
[Pxy, f_pxy] = cpsd(x, y, window, overlap, nfft, fs);
phase_xy = angle(Pxy);

% Phase relationship
subplot(2,1,1);
freq_range_phase = (f_pxy >= 0.1) & (f_pxy <= 8);
plot(f_pxy(freq_range_phase), phase_xy(freq_range_phase)*180/pi, 'r-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Phase Difference (degrees)');
title('Phase Relationship between Signals');
grid on;
xlim([0.1, 8]);

% Cross-spectral density magnitude
subplot(2,1,2);
plot(f_pxy(freq_range_phase), abs(Pxy(freq_range_phase)), 'k-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Cross-Spectral Density Magnitude');
title('Cross-Spectral Density');
grid on;
xlim([0.1, 8]);

%% Summary statistics
fprintf('\n=== Coherence Analysis Summary ===\n');
fprintf('Sampling rate: %d Hz\n', fs);
fprintf('Data length: %.1f minutes (%d samples)\n', duration/60, N);
fprintf('Frequency resolution: %.4f Hz\n', fs/window_length);
fprintf('Window length: %.1f seconds\n', window_length/fs);
fprintf('Number of segments: %d\n', num_segments);
fprintf('\nMean coherence by frequency band:\n');
for i = 1:length(band_names)
    fprintf('%s Hz: %.3f\n', ...
        sprintf('%.1f-%.1f', bands.(band_names{i})(1), bands.(band_names{i})(2)), ...
        band_coherence(i));
end
confidence_level = empirical_threshold;
% Find peaks in coherence
[peaks, peak_locs] = findpeaks(Cxy_roi, f_roi, 'MinPeakHeight', confidence_level);
if ~isempty(peaks)
    fprintf('\nSignificant coherence peaks:\n');
    for i = 1:length(peaks)
        fprintf('%.3f Hz: Coherence = %.3f\n', peak_locs(i), peaks(i));
    end
else
    fprintf('\nNo significant coherence peaks found above confidence threshold.\n');
end
