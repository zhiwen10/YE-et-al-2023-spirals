%%
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'YE-et-al-2023-spirals')));            % paper repository
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'npy-matlab'))) % kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\wheelAnalysis'));
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);                                             % atlas brain boundry binary mask 
[row,col] = find(BW);
brain_index = [col,row];                                                   % position idenx within the brain boundry
%% load session table
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
T = T(T.face &T.eye_DLC,:);
%%
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
%% coherence amp and phase
% behavior = 'face';
behavior = 'pupil';
[f_pxy,Cxy_all,Pxy_all] = coherence2traces(T,data_folder,pixel,behavior);
% surrogate threshold
load([behavior '_coherence_data_7sessions.mat'],'pooled_surrogate_threshold',...
    'freq_axis');
n_sessions = 7;
freq_mask2 = (freq_axis >= 0.05) & (freq_axis <= 8);
pooled_surrogate_threshold_roi = pooled_surrogate_threshold(freq_mask2);
%%
Cxy_mean = mean(Cxy_all,2);
Cxy_sem = std(Cxy_all,[],2)./sqrt(size(Cxy_all,2));
min_coh = Cxy_mean -Cxy_sem;
min_coh = min_coh(freq_mask2);
a = (min_coh>pooled_surrogate_threshold_roi);
index = find(a==0,1,'first');
freq_roi = freq_axis(freq_mask2);
freq1 = freq_roi(index);
% phase mean and sem
index2 = find(f_pxy == freq1);
freq_mask3 = freq_mask2;
freq_mask3(index2+1:end) = 0;
phase_xy_all = angle(Pxy_all);
phase_xy_roi2 = phase_xy_all(freq_mask3,:);
phase_xy_roi2_mean_freq_range = mean(phase_xy_roi2,1);
phase_xy_roi2_mean = mean(phase_xy_roi2_mean_freq_range)./pi*180;
phase_xy_roi2_sem = std(phase_xy_roi2_mean_freq_range)./sqrt(size(phase_xy_roi2,2))./pi*180;
[hh1,pp1] = ttest(phase_xy_roi2_mean_freq_range,0);
%% phase-amplitude coupling
n_bins = 18; 
low_freq_band = [0.05, 0.5];                                               % Phase-providing frequency (slow oscillation)
high_freq_band = [2, 8];                                                   % Amplitude-modulated frequency (fast oscillation)
[phase_centers,mean_amp_per_bin_all2] = phaseAmpCoupling(T,data_folder,...
    pixel,behavior, n_bins,low_freq_band,high_freq_band);
mean_amp_per_bin_all2 = mean_amp_per_bin_all2*100;                         % df/f convert to percentage
amp_per_bin_mean2 = mean(mean_amp_per_bin_all2,3);                         % average across all pixels
% significance test (circular-linear test)
phase_centers2 = repmat(phase_centers,7,1)';
[rho pval] = circ_corrcl(phase_centers2(:), amp_per_bin_mean2(:));
% ttest
amp_high = amp_per_bin_mean2(1,:);
amp_high_mean = mean(amp_high); amp_high_sem = std(amp_high)./sqrt(numel(amp_high));
amp_low = amp_per_bin_mean2(10,:);
amp_low_mean = mean(amp_low); amp_low_sem = std(amp_low)./sqrt(numel(amp_low));
[ha1,pa1] = ttest(amp_high,amp_low);
%% phase-spiral histogram
spirals_sort = phaseSpiralHistogram(T,data_folder,behavior,...
    n_bins,brain_index,low_freq_band);
count_sample = spiralDensityBins(T,data_folder,spirals_sort);
%%
count_sample2 = count_sample';
[rho2 pval2] = circ_corrcl(phase_centers2(:), count_sample2(:));
% ttest
count_high = count_sample2(1,:);
count_high_mean = mean(count_high); count_high_sem = std(count_high)./sqrt(numel(count_high));
count_low = count_sample2(10,:);
count_low_mean = mean(count_low); count_low_sem = std(count_low)./sqrt(numel(count_low));
[ha2,pa2] = ttest(count_high,count_low);
%% Plotting results
h1 = figure('Position', [100, 100, 800, 600]);
freq_range_phase = (f_pxy >= 0.05) & (f_pxy <= 8);
color1 = cbrewer2('seq','Greys',10);
color2 = '-b';
ax1 = subplot(2,2,1);
for session = 1:size(T,1)
    Cxy = Cxy_all(:,session);
    plot(log10(f_pxy(freq_range_phase)), Cxy(freq_range_phase),'color',color1(session+3,:), 'LineWidth', 1);
    xlabel('Frequency (Hz) - Log Scale');
    ylabel('Coherence');
    title('Coherence - Log Frequency Scale');
    xlim([log10(0.05), log10(8)]);
    ylim([0, 1]);
    xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
    xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
    hold on;
end
shadedErrorBar(log10(f_pxy(freq_range_phase)), Cxy_mean(freq_range_phase),...
    Cxy_sem(freq_range_phase), 'lineprops', color2);
hold on;
plot([log10(0.05),log10(freq1)],[0.74,0.74],'k','linewidth',4);
hold on;
plot(log10(freq_roi), pooled_surrogate_threshold(freq_mask2), 'r--', 'LineWidth', 2);
legend('Mean ± SEM', 'Surrogate 95%', 'Location', 'best');
%
ax2 = subplot(2,2,2);
freq_range_phase = (f_pxy >= 0.05) & (f_pxy <= 8);
for session = 1:size(T,1)
    Cxy = Cxy_all(:,session);
    Pxy = Pxy_all(:,session);
    phase_xy = angle(Pxy);
    % Phase relationship
    plot(log10(f_pxy(freq_range_phase)), phase_xy(freq_range_phase),'color',color1(session+3,:), 'LineWidth', 1);
    hold on;
end
phase_xy_all = angle(Pxy_all);
phase_mean = mean(phase_xy_all,2);
phase_sem = std(phase_xy_all,[],2)./sqrt(size(phase_xy_all,2));
shadedErrorBar(log10(f_pxy(freq_range_phase)), phase_mean(freq_range_phase),...
    phase_sem(freq_range_phase), 'lineprops', color2);
yline(0,'--r');
xlim([log10(0.05), log10(8)]);
xticks(log10([0.05,0.1,0.2,0.5,1,2,4,8]));
xticklabels([0.05,0.1,0.2,0.5,1,2,4,8]);
yticks([-pi:pi/2:pi]);
yticklabels({'-pi','-pi/2','0','pi/2','pi'});
xlabel('Frequency (Hz) - Log Scale');
ylabel('Phase Difference');
title('Phase Relationship between Signals');

ax3 = subplot(2,2,3);
amp_per_bin_mean3 = mean(amp_per_bin_mean2,2);
amp_per_bin_sem3 = std(amp_per_bin_mean2,[],2)./sqrt(size(amp_per_bin_mean2,2));
for session = 1:size(T,1)
    plot(phase_centers, amp_per_bin_mean2(:,session),'color',color1(session+3,:));
    hold on;
end
shadedErrorBar(phase_centers, amp_per_bin_mean3,amp_per_bin_sem3, 'lineprops', color2)
xlim([-pi,pi]);
ylim([-0.2,0.2]);
yticks([-0.2:0.1:0.2]);
xticks([-pi:pi/2:pi]);
xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
xlabel('Phase of behavior (0.1-0.5 Hz)');
ylabel('2-8 Hz amp (dF/F,%)');
title('Phase-amplitude coulping');

ax4 = subplot(2,2,4);
for session = 1:size(count_sample,1)
    spiral = squeeze(count_sample(session,:));     
    plot(phase_centers,spiral,'color',color1(session+3,:));
    hold on;
end
xticks([-pi:pi/2:pi]);
xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
xlim([-pi,pi]);
ylim([0,6]);
yticks([0:2:6]);
count_mean = mean(count_sample,1);
count_sem = std(count_sample,[],1)./sqrt(size(count_sample,1));
hold on;
shadedErrorBar(phase_centers, count_mean, count_sem, 'lineprops', color2);
xlabel('Phase of behavior (0.1-0.5 Hz)');
ylabel('Peak density (spirals/mm2*s)');
title('Phase-rotating wave density coulping');

print(h1, 'pupil_summary2.pdf','-dpdf', '-bestfit', '-painters');
%%
h2 = figure('Position', [100, 100, 800, 200]);
color3 = colorcet('C4','N',181);
color3 = circshift(color3,45,1);
ax1 = subplot(1,1,1);
x = -pi:pi/90:pi;
y = cos(x);
plot(x,y);
hold on;
x1 = x(1:end-1);
x2 = x(2:end);
y1 = y(1:end-1);
y2 = y(2:end);
for j = 1:numel(x)-1
    plot([x1(j),x2(j)],[y1(j),y2(j)],'linewidth',4,'color',color3(j,:));
    hold on;
end
yline(0);
xlim([-pi,pi]);
print(h2, 'arousal_cycle.pdf','-dpdf', '-bestfit', '-painters');