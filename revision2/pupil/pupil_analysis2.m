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
%     tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
%     pupil_mean = readNPY(tFile1); 
%     % pupilsize = [tUp,pupil_mean];
%     pupil_mean2 = pupil_mean(1:2:end);
%     if numel(pupil_mean2)<size(t,1)
%         pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
%     elseif numel(pupil_mean2)>size(t,1)
%         pupil_mean2 = pupil_mean2(1:numel(t));
%     end 
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
    %%
    freq_range = (f >= 0.05) & (f <= 8);
    f_roi = f(freq_range);
    Cxy_roi = Cxy(freq_range);
    % Compute cross-spectral density and phase
    [Pxy, f_pxy] = cpsd(x, y, window, overlap, nfft, fs);
    phase_xy = angle(Pxy);
    %%
    f_pxy_all(:,kk) = f_pxy;
    Cxy_all(:,kk) = Cxy;
    Pxy_all(:,kk) = Pxy;
end
%% Plotting results
h1 = figure('Position', [100, 100, 900, 300]);
freq_range_phase = (f_pxy >= 0.05) & (f_pxy <= 8);
for kk = 1:size(T,1)
    f_pxy = f_pxy_all(:,kk);
    Cxy = Cxy_all(:,kk);
    Pxy = Pxy_all(:,kk);
    phase_xy = angle(Pxy);
    subplot(1,2,1);
    plot(log10(f_pxy(freq_range_phase)), Cxy(freq_range_phase), 'k-', 'LineWidth', 1, 'DisplayName', 'Observed');
    xlabel('Frequency (Hz) - Log Scale');
    ylabel('Coherence');
    title('Coherence - Log Frequency Scale');
    xlim([log10(0.05), log10(8)]);
    ylim([0, 1]);
    xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
    xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
    hold on;
    % Additional analysis: Phase and cross-spectral density

    % Phase relationship
    subplot(1,2,2);
    plot(log10(f_pxy(freq_range_phase)), phase_xy(freq_range_phase), 'k-', 'LineWidth', 1);
    xlabel('Frequency (Hz) - Log Scale');
    ylabel('Phase Difference');
    title('Phase Relationship between Signals');
    xlim([log10(0.05), log10(8)]);
    xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
    xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
    hold on;
end

Cxy_mean = mean(Cxy_all,2);
Cxy_sem = std(Cxy_all,[],2)./sqrt(size(Cxy_all,2));
subplot(1,2,1);
errorbar(log10(f_pxy(freq_range_phase)), Cxy_mean(freq_range_phase),...
    Cxy_sem(freq_range_phase),'b-', 'LineWidth', 1, 'DisplayName', 'Observed');
subplot(1,2,2);
phase_xy_all = angle(Pxy_all);
phase_mean = mean(phase_xy_all,2);
phase_sem = std(phase_xy_all,[],2)./sqrt(size(phase_xy_all,2));
errorbar(log10(f_pxy(freq_range_phase)), phase_mean(freq_range_phase),...
    phase_sem(freq_range_phase),'b-', 'LineWidth', 1, 'DisplayName', 'Observed');
yline(0,'--r');
yticks([-pi:pi/2:pi]);
yticklabels({'-pi','-pi/2','0','pi/2','pi'});
%%
h1= figure('Renderer', 'painters', 'Position', [100 100 400 400]);
subplot(1,1,1);
freq_range_phase = (f_pxy >= 0.05) & (f_pxy <= 8);
for kk = 1:size(T,1)
    f_pxy = f_pxy_all(:,kk);
    Cxy = Cxy_all(:,kk);
    Pxy = Pxy_all(:,kk);
    phase_xy = angle(Pxy);
    % Phase relationship
    plot(log10(f_pxy(freq_range_phase)), phase_xy(freq_range_phase), 'k-', 'LineWidth', 1);
    xlabel('Frequency (Hz) - Log Scale');
    ylabel('Phase Difference');
    title('Phase Relationship between Signals');
    xlim([log10(0.05), log10(8)]);
    xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
    xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
    hold on;
end

Cxy_mean = mean(Cxy_all,2);
Cxy_sem = std(Cxy_all,[],2)./sqrt(size(Cxy_all,2));
phase_xy_all = angle(Pxy_all);
phase_mean = mean(phase_xy_all,2);
phase_sem = std(phase_xy_all,[],2)./sqrt(size(phase_xy_all,2));
% errorbar(log10(f_pxy(freq_range_phase)), phase_mean(freq_range_phase),...
%     phase_sem(freq_range_phase),'b-', 'LineWidth', 1);
shadedErrorBar(log10(f_pxy(freq_range_phase)), phase_mean(freq_range_phase),...
    phase_sem(freq_range_phase), 'lineprops', '-b');
yline(0,'--r');
yticks([-pi:pi/2:pi]);
yticklabels({'-pi','-pi/2','0','pi/2','pi'});
print(h1, 'motion_phase_diff_7sessions.pdf',...
    '-dpdf', '-bestfit', '-painters');