%% Multi-Session Coherence Significance Testing
% Determining significant frequencies across 15 sessions using surrogate methods
%%
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
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
% T = T(T.face &T.eye_DLC,:);
T = T(logical(T.face),:);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% Parameters
fs = 35;                    % Sampling frequency
n_sessions = 13;            % Number of sessions
% duration = 30 * 60;         % 30 minutes
% N = fs * duration;          % Samples per session

% Coherence parameters
window_length = 2^12;       % 4096 samples for good freq resolution
overlap = window_length/2;   % 50% overlap
nfft = window_length;
window = hamming(window_length);

% Surrogate parameters
n_surrogates = 500;         % Per session (computational balance)
freq_range_interest = [0.1, 8]; % Focus frequency range

% Multiple comparison correction methods
alpha = 0.05;

%% Initialize storage
coherences_all = cell(n_sessions, 1);
surr_coherences_all = cell(n_sessions, 1);
freq_axis = [];

fprintf('Processing %d sessions for coherence significance testing...\n', n_sessions);

%% Session-by-session analysis
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
    trace = squeeze(Utransformed(pixel(7,1),pixel(7,2),1:50))'*V(1:50,:);
    %%
    load(fullfile(data_folder,'spirals','spirals_index',...
        [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    %% Parameters
    fs = 35;                    % Sampling frequency (Hz)
    duration = numel(t)/35;         % 30 minutes in seconds
    N = fs * duration;          % Total number of samples (63,000)
    % Assume your data is loaded as:
    % x = pupil_mean2;   % Replace with your actual data
    x = image_energy2;
    y = double(trace)';  % Replace with your actual data
    % y = pupil_mean2;
    sizeN = min(numel(x),numel(y));
    x = x(1:sizeN);
    y = y(1:sizeN);
    %%
    % Ensure column vectors
    if size(x, 2) > size(x, 1), x = x'; end
    if size(y, 2) > size(y, 1), y = y'; end
    
    %% Compute observed coherence
    [Cxy, f] = mscohere(x, y, window, overlap, nfft, fs);
    coherences_all{kk} = Cxy;
    
    % Store frequency axis (same for all sessions)
    if kk == 1
        freq_axis = f;
        n_freqs = length(f);
        
        % Pre-allocate arrays now that we know frequency vector size
        observed_coherences = zeros(n_sessions, n_freqs);
        surrogate_thresholds = zeros(n_sessions, n_freqs);
        session_pvalues = zeros(n_sessions, n_freqs);
        significant_freqs_per_session = false(n_sessions, n_freqs);
    end
    
    observed_coherences(kk, :) = Cxy;
    
    %% Generate surrogate data for this session
    surrogate_coherences_session = zeros(n_surrogates, n_freqs);
    
    for surr = 1:n_surrogates
        % Phase randomization of x signal
        X_fft = fft(x);
        random_phases = 2*pi*rand(size(X_fft));
%         % Ensure conjugate symmetry for real output
%         random_phases(1) = 0; % DC component
%         if mod(length(x), 2) == 0
%             random_phases(end/2+1) = 0; % Nyquist frequency for even length
%         end
%         % Make phases conjugate symmetric
%         random_phases(ceil(end/2)+1:end) = -random_phases(ceil(end/2):-1:2);
%         
        x_surrogate = real(ifft(abs(X_fft) .* exp(1i*random_phases)));
        
        % Compute surrogate coherence
        [Cxy_surr, ~] = mscohere(x_surrogate, y, window, overlap, nfft, fs);
        surrogate_coherences_session(surr, :) = Cxy_surr;
    end
    
    % Store surrogate data
    surr_coherences_all{kk} = surrogate_coherences_session;
    
    %% Significance testing for this session
    % 95th percentile threshold for each frequency
    surrogate_thresholds(kk, :) = prctile(surrogate_coherences_session, 95, 1);
    
    % P-values: proportion of surrogates >= observed
    for freq = 1:n_freqs
        session_pvalues(kk, freq) = sum(surrogate_coherences_session(:, freq) >= Cxy(freq)) / n_surrogates;
    end
    
    % Identify significant frequencies (uncorrected)
    significant_freqs_per_session(kk, :) = session_pvalues(kk, :) < alpha;
end

%% Group-level significance testing across sessions

% Method 1: Fisher's method for combining p-values across sessions
fisher_chi2 = -2 * sum(log(session_pvalues + eps), 1); % Add eps to avoid log(0)
fisher_pvalues = 1 - chi2cdf(fisher_chi2, 2*n_sessions);

% Method 2: Proportion of significant sessions per frequency
proportion_significant = sum(significant_freqs_per_session, 1) / n_sessions;

% Method 3: Mean coherence across sessions vs pooled surrogate distribution
mean_coherence_across_sessions = mean(observed_coherences, 1);
std_coherence_across_sessions = std(observed_coherences, 0, 1);

% Create pooled surrogate distribution
pooled_surrogates = [];
for kk = 1:n_sessions
    pooled_surrogates = [pooled_surrogates; surr_coherences_all{kk}];
end
pooled_surrogate_threshold = prctile(pooled_surrogates, 95, 1);

% Method 4: Binomial test for each frequency
% Test if proportion of significant sessions > chance (5%)
binomial_pvalues = zeros(1, n_freqs);
for freq = 1:n_freqs
    n_sig = sum(significant_freqs_per_session(:, freq));
    binomial_pvalues(freq) = 1 - binocdf(n_sig-1, n_sessions, alpha);
end

%% Multiple comparison corrections
% Focus on frequency range of interest
freq_mask = (freq_axis >= freq_range_interest(1)) & (freq_axis <= freq_range_interest(2));
n_freqs_tested = sum(freq_mask);

% Bonferroni correction
fisher_pvalues_bonf = fisher_pvalues * n_freqs_tested;
fisher_pvalues_bonf(fisher_pvalues_bonf > 1) = 1;

% False Discovery Rate (FDR) correction using Benjamini-Hochberg
fisher_pvalues_fdr = mafdr(fisher_pvalues(freq_mask), 'BHFDR', true);
fisher_pvalues_fdr_full = ones(size(fisher_pvalues));
fisher_pvalues_fdr_full(freq_mask) = fisher_pvalues_fdr;

% Final significance calls
sig_fisher_uncorrected = fisher_pvalues < alpha;
sig_fisher_bonferroni = fisher_pvalues_bonf < alpha;
sig_fisher_fdr = fisher_pvalues_fdr_full < alpha;
sig_binomial = binomial_pvalues < alpha;
sig_proportion = proportion_significant > 0.5; % Majority of sessions

%% Visualization
figure('Position', [50, 50, 1800, 1200]);

% Plot 1: Individual session coherences
subplot(3, 4, 1);
freq_roi = freq_axis(freq_mask);
coherences_roi = observed_coherences(:, freq_mask);
plot(freq_roi, coherences_roi', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold on;
plot(freq_roi, mean(coherences_roi, 1), 'b-', 'LineWidth', 3);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('All Sessions + Mean');
grid on;
xlim(freq_range_interest);

% Plot 2: Mean coherence with error bars
subplot(3, 4, 2);
mean_coh_roi = mean(coherences_roi, 1);
sem_coh_roi = std(coherences_roi, 0, 1) / sqrt(n_sessions);
errorbar(freq_roi, mean_coh_roi, sem_coh_roi, 'b-', 'LineWidth', 2);
hold on;
plot(freq_roi, pooled_surrogate_threshold(freq_mask), 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Coherence ± SEM');
title('Group Mean vs Surrogate Threshold');
legend('Mean ± SEM', 'Surrogate 95%', 'Location', 'best');
grid on;
xlim(freq_range_interest);

% Plot 3: Proportion of significant sessions
subplot(3, 4, 3);
bar(freq_roi, proportion_significant(freq_mask), 'FaceColor', [0.6 0.8 0.6]);
hold on;
yline(0.5, 'r--', 'LineWidth', 2, 'DisplayName', '50% threshold');
xlabel('Frequency (Hz)');
ylabel('Proportion Significant');
title('Fraction of Sessions Significant');
legend('Location', 'best');
grid on;
xlim(freq_range_interest);
ylim([0, 1]);

% Plot 4: Fisher's combined p-values
subplot(3, 4, 4);
semilogy(freq_roi, fisher_pvalues(freq_mask), 'b-', 'LineWidth', 2);
hold on;
semilogy(freq_roi, fisher_pvalues_fdr_full(freq_mask), 'g-', 'LineWidth', 2);
yline(alpha, 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('P-value (log scale)');
title('Fisher Combined P-values');
legend('Uncorrected', 'FDR corrected', 'α = 0.05', 'Location', 'best');
grid on;
xlim(freq_range_interest);

% Plot 5-7: Heatmaps for different significance methods
subplot(3, 4, 5);
imagesc(freq_roi, 1:n_sessions, significant_freqs_per_session(:, freq_mask));
colormap(gca, [1 1 1; 1 0 0]);
xlabel('Frequency (Hz)');
ylabel('Session');
title('Individual Session Significance');
colorbar('Ticks', [0, 1], 'TickLabels', {'Not Sig', 'Significant'});

subplot(3, 4, 6);
sig_map_fisher = double(sig_fisher_fdr(freq_mask));
imagesc(freq_roi, 1, sig_map_fisher);
colormap(gca, [1 1 1; 0 1 0]);
xlabel('Frequency (Hz)');
title('Fisher FDR Significant');
set(gca, 'YTick', []);

subplot(3, 4, 7);
sig_map_proportion = double(sig_proportion(freq_mask));
imagesc(freq_roi, 1, sig_map_proportion);
colormap(gca, [1 1 1; 0 0 1]);
xlabel('Frequency (Hz)');
title('Majority Sessions Significant');
set(gca, 'YTick', []);

% Plot 8: Summary bar chart of significant frequencies
subplot(3, 4, 8);
methods = {'Individual', 'Fisher FDR', 'Proportion>50%', 'Binomial'};
n_sig_freqs = [sum(any(significant_freqs_per_session(:, freq_mask), 1)), ...
               sum(sig_fisher_fdr(freq_mask)), ...
               sum(sig_proportion(freq_mask)), ...
               sum(sig_binomial(freq_mask))];

bar(n_sig_freqs, 'FaceColor', [0.8 0.6 0.4]);
set(gca, 'XTickLabel', methods, 'XTickLabelRotation', 45);
ylabel('Number of Significant Frequencies');
title(sprintf('Significant Frequencies\n(out of %d tested)', n_freqs_tested));
grid on;

% Plots 9-12: Individual session examples
for i = 1:4
    subplot(3, 4, 8+i);
    sess_idx = i * 3 + 1; % Show sessions 4, 7, 10, 13
    if sess_idx <= n_sessions
        plot(freq_roi, coherences_roi(sess_idx, :), 'b-', 'LineWidth', 2);
        hold on;
        plot(freq_roi, surrogate_thresholds(sess_idx, freq_mask), 'r--', 'LineWidth', 1.5);
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        title(sprintf('Session %d Example', sess_idx));
        grid on;
        xlim(freq_range_interest);
        ylim([0, 1]);
    end
end

sgtitle('Multi-Session Coherence Significance Analysis (15 Sessions)', 'FontSize', 16);
%%
figure;
% Plot 2: Mean coherence with error bars
subplot(1,1,1);
mean_coh_roi = mean(coherences_roi, 1);
sem_coh_roi = std(coherences_roi, 0, 1) / sqrt(n_sessions);
errorbar(log10(freq_roi), mean_coh_roi, 3*sem_coh_roi, 'b-', 'LineWidth', 2);
hold on;
plot(log10(freq_roi), pooled_surrogate_threshold(freq_mask), 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Coherence ± SEM');
title('Group Mean vs Surrogate Threshold');
legend('Mean ± SEM', 'Surrogate 95%', 'Location', 'best');
grid on;
% xlim(freq_range_interest);
xticks(log10([0.1,0.2,0.5,1,2,4,8]));
xticklabels([0.1,0.2,0.5,1,2,4,8]);
%%
figure;
% significant_freqs_per_session2 = double(significant_freqs_per_session);
subplot(1,1,1);
for i = 1:13
    scatter(log10(freq_roi), significant_freqs_per_session(i, freq_mask)*i,'|k');
    hold on;
%     xlabel('Frequency (Hz)');
%     ylabel('Session');
%     title('Individual Session Significance');
%     colorbar('Ticks', [0, 1], 'TickLabels', {'Not Sig', 'Significant'});
end
xticks(log10([0.1,0.2,0.5,1,2,4,8]));
xticklabels([0.1,0.2,0.5,1,2,4,8]);
%%
figure;
% Plot 2: Mean coherence with error bars
subplot(1,1,1);
freq_roi = freq_axis(freq_mask);
coherences_roi = observed_coherences(:, freq_mask);
plot(log10(freq_roi), coherences_roi', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold on;
mean_coh_roi = mean(coherences_roi, 1);
sem_coh_roi = std(coherences_roi, 0, 1) / sqrt(n_sessions);
shadedErrorBar(log10(freq_roi), mean_coh_roi, 3*sem_coh_roi, 'lineprops', '-g');
hold on;
plot(log10(freq_roi), pooled_surrogate_threshold(freq_mask), 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz) - Log Scale');
ylabel('Coherence');
xticks(log10([0.1,0.2,0.5,1,2,4,8]));
xticklabels([0.1,0.2,0.5,1,2,4,8]);
box off;
pooled_surrogate_threshold_roi = pooled_surrogate_threshold(freq_mask);
min_coh = mean_coh_roi-3*sem_coh_roi;
a = (min_coh>pooled_surrogate_threshold_roi);
index = find(a==0,1,'first');
freq1 = freq_roi(index);
hold on;
plot([log10(0.1),log10(freq1)],[0.74,0.74],'k','linewidth',4);
legend('Mean ± 3*SEM', 'Surrogate 95%', 'Location', 'best');