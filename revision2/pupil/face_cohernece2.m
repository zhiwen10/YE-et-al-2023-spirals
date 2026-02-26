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
T = T(T.face &T.eye_DLC,:);
% T = T(logical(T.face),:);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% Parameters
fs = 35;                    % Sampling frequency
n_sessions = 7;            % Number of sessions
% duration = 30 * 60;         % 30 minutes
% N = fs * duration;          % Samples per session
% Coherence parameters
window_length = 2^12;       % 4096 samples for good freq resolution
overlap = window_length/2;   % 50% overlap
nfft = window_length;
window = hamming(window_length);
% Surrogate parameters
n_surrogates = 500;         % Per session (computational balance)
freq_range_interest = [0.01, 8]; % Focus frequency range
% Multiple comparison correction methods
alpha = 0.05;
%% Initialize storage
coherences_all = cell(n_sessions, 1);
surr_coherences_all = cell(n_sessions, 1);
freq_axis = [];
fprintf('Processing %d sessions for coherence significance testing...\n', n_sessions);
% Session-by-session analysis
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
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    trace = squeeze(Utransformed(pixel(7,1),pixel(7,2),1:50))'*V(1:50,:);
    %%
%     load(fullfile(data_folder,'spirals','spirals_index',...
%         [fname '_motion_energy.mat']));
%     image_energy2(isnan(image_energy2)) = 0;
    %%
    tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
    pupil_mean = readNPY(tFile1); 
    pupil_mean2 = pupil_mean(1:2:end);
    if numel(pupil_mean2)<size(t,1)
        pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
    elseif numel(pupil_mean2)>size(t,1)
        pupil_mean2 = pupil_mean2(1:numel(t));
    end 
    %% Parameters
    fs = 35;                    % Sampling frequency (Hz)
    duration = numel(t)/35;         % 30 minutes in seconds
    N = fs * duration;          % Total number of samples (63,000)
    % Assume your data is loaded as:
    x = pupil_mean2;   % Replace with your actual data
    % x = image_energy2;
    y = double(trace)';  % Replace with your actual data
    sizeN = min(numel(x),numel(y));
    x = x(1:sizeN);
    y = y(1:sizeN);
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
% Create pooled surrogate distribution
pooled_surrogates = [];
for kk = 1:n_sessions
    pooled_surrogates = [pooled_surrogates; surr_coherences_all{kk}];
end
pooled_surrogate_threshold = prctile(pooled_surrogates, 95, 1);
%%
freq_mask = (freq_axis >= freq_range_interest(1)) & (freq_axis <= freq_range_interest(2));
%%
save('pupil_coherence_data.mat','pooled_surrogate_threshold',...
    'session_pvalues','significant_freqs_per_session','surr_coherences_all',...
    'surrogate_thresholds','freq_mask','freq_axis','observed_coherences',...
    'n_sessions','pooled_surrogate_threshold');
%%
load('motion_coherence_data2.mat','pooled_surrogate_threshold',...
    'session_pvalues','significant_freqs_per_session','surr_coherences_all',...
    'surrogate_thresholds','freq_mask','freq_axis','observed_coherences',...
    'n_sessions','pooled_surrogate_threshold');

%%
n_sessions = 7;
freq_mask2 = (freq_axis >= 0.05) & (freq_axis <= 8);
h1= figure('Renderer', 'painters', 'Position', [100 100 400 400]);
% Plot 2: Mean coherence with error bars
subplot(1,1,1);
freq_roi = freq_axis(freq_mask2);
coherences_roi = observed_coherences(1:n_sessions, freq_mask2);
plot(log10(freq_roi), coherences_roi', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold on;
mean_coh_roi = mean(coherences_roi, 1);
sem_coh_roi = std(coherences_roi, 0, 1) / sqrt(n_sessions);
shadedErrorBar(log10(freq_roi), mean_coh_roi, sem_coh_roi, 'lineprops', '-g');
hold on;
plot(log10(freq_roi), pooled_surrogate_threshold(freq_mask2), 'r--', 'LineWidth', 2);
xlabel('Frequency (Hz) - Log Scale');
ylabel('Coherence');
xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
xlim([log10(0.05),log10(8)]);
box off;
pooled_surrogate_threshold_roi = pooled_surrogate_threshold(freq_mask2);
min_coh = mean_coh_roi-sem_coh_roi;
a = (min_coh>pooled_surrogate_threshold_roi);
index = find(a==0,1,'first');
freq1 = freq_roi(index);
hold on;
plot([log10(0.05),log10(freq1)],[0.74,0.74],'k','linewidth',4);
legend('Mean ± SEM', 'Surrogate 95%', 'Location', 'best');
print(h1, 'motion_coherence_7sessions.pdf',...
    '-dpdf', '-bestfit', '-painters');