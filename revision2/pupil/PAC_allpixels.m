%%
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\wheelAnalysis'));
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
T = T(logical(T.eye_DLC),:);
% T = T(logical(T.face),:);
% T = T([1,2,6:13],:);
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
low_freq_band = [0.05, 0.5]; % Phase-providing frequency (slow oscillation)
high_freq_band = [2, 8];    % Amplitude-modulated frequency (fast oscillation)
% high_freq_band = [0.5, 2];    % Amplitude-modulated frequency (fast oscillation)
mean_amp_per_bin_all = [];
mean_amp_per_bin_all2 = [];
mean_amp_per_bin_all_sum = [];
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
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    %%
    tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
    pupil_mean = readNPY(tFile1); 
    pupil_mean2 = pupil_mean(1:2:end);
    if numel(pupil_mean2)<size(t,1)
        pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
    elseif numel(pupil_mean2)>size(t,1)
        pupil_mean2 = pupil_mean2(1:numel(t));
    end 
    signal2 = pupil_mean2;
    %% rotaryEncoder
%     sigName = 'rotaryEncoder';
%     tlFile = fullfile(serverRoot, [sigName '.raw.npy']); 
%     pd = readNPY(tlFile);
%     fs1 = 35;
%     wh = wheel.correctCounterDiscont(pd);
%     tlFile = fullfile(serverRoot, [sigName '.timestamps_Timeline.npy']);
%     tlTimes = readNPY(tlFile);
%     tt = tsToT(tlTimes, numel(pd));  
%     fs1 = 1/mean(diff(tt));
%     wh = wheel.correctCounterDiscont(pd);
%     vel = wheel.computeVelocity2(wh, 0.01, fs1); 
%     vel2 = interp1(tt,vel,t);
%     signal2 = vel2;
    %%
%     load(fullfile(data_folder,'spirals','spirals_index',...
%     [fname '_motion_energy.mat']));
%     image_energy2(isnan(image_energy2)) = 0;
%     if numel(image_energy2)<size(t,1)
%         image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
%     elseif numel(image_energy2)>size(t,1)
%         image_energy2 = image_energy2(1:numel(t));
%     end 
%     signal2 = image_energy2;
    %% Method 1: Modulation Index (Tort et al., 2010) - Most common approach
    % Design filters for frequency bands
    % Low frequency (phase) - 4th order Butterworth
    [b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
    signal_low = filtfilt(b_low, a_low, signal2);    
    % Extract phase from low frequency signal using Hilbert transform
    phase_low = angle(hilbert(signal_low));
    %%
    [b_high, a_high] = butter(4, high_freq_band/(fs/2), 'bandpass');
    for j = 1:7
        %%
        trace = squeeze(Utransformed(pixel(j,1),pixel(j,2),1:50))'*V(1:50,:);
        signal = double(trace)'./mimgtransformed(pixel(j,1),pixel(j,2));
        %%
%         trace3 = squeeze(Utransformed(pixel(7,1),pixel(7,2),1:50))'*V(1:50,:);
%         signal3 = double(trace3)'./mimgtransformed(pixel(7,1),pixel(7,2));
%         figure;
%         plot(t,signal2/1000000,'r');
%         hold on;
%         plot(t,vel2/10000,'g');
%         hold on;
%         plot(t,signal,'k');
%         hold on;
%         plot(t,signal3-0.1,'k');
        %%
        % High frequency (amplitude) - 4th order Butterworth  
        signal_high = filtfilt(b_high, a_high, signal);
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
        MI(:,kk,j) = (H_max - H) / H_max; % Modulation Index (0 = no coupling, 1 = perfect coupling)
        mean_amp_per_bin_all(:,kk,j) = mean_amp_per_bin;
    end
end
%%
mean_amp_per_bin = mean(mean_amp_per_bin_all,1);
mean_amp_per_bin_all2 = mean_amp_per_bin_all-mean_amp_per_bin;
mean_amp_per_bin_all2 = mean_amp_per_bin_all2*100;
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 800 600]);
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
for j = 1:7
    ax1 = subplot(2,4,j);
    amp_per_bin1 = squeeze(mean_amp_per_bin_all2(:,:,j));
    amp_per_bin_mean = mean(amp_per_bin1,2);
    amp_per_bin_sem = std(amp_per_bin1,[],2)./sqrt(size(amp_per_bin1,2));
    plot(phase_centers, amp_per_bin1,'k');
%     hold on;
%     errorbar(phase_centers, amp_per_bin_mean,amp_per_bin_sem,'r');
    hold on;
    shadedErrorBar(phase_centers, amp_per_bin_mean,amp_per_bin_sem, 'lineprops', '-r')
    xlim([-pi,pi]);
    ylim([-0.8,0.8]);
    xticks([-pi:pi/2:pi]);
    xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
    xlabel('Phase of facial motion (0.1-0.5 Hz)');
    ylabel('Mean High Freq Amplitude  (dF/F,%)');
end
print(h1, 'pupil_phase_amp_coupling_sessions.pdf','-dpdf', '-bestfit', '-painters');
%%
h2 = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
amp_per_bin_mean2 = mean(mean_amp_per_bin_all2,3);
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
ax1 = subplot(1,1,1);
amp_per_bin_mean3 = mean(amp_per_bin_mean2,2);
amp_per_bin_sem3 = std(amp_per_bin_mean2,[],2)./sqrt(size(amp_per_bin_mean2,2));
plot(phase_centers, amp_per_bin_mean2,'k');
% hold on;
% errorbar(phase_centers, amp_per_bin_mean3,amp_per_bin_sem3,'r');
hold on;
shadedErrorBar(phase_centers, amp_per_bin_mean3,amp_per_bin_sem3, 'lineprops', '-r')
xlim([-pi,pi]);
ylim([-0.2,0.2]);
yticks([-0.2:0.1:0.2]);
xticks([-pi:pi/2:pi]);
xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
xlabel('Phase of facial motion (0.1-0.5 Hz)');
ylabel('Mean High Freq Amplitude (dF/F,%)');
print(h2, 'pupil_phase_amp_coupling_all.pdf','-dpdf', '-bestfit', '-painters');
%%
ratio = amp_per_bin_mean2(4,:)-amp_per_bin_mean2(13,:);