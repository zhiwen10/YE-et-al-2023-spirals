%% behavior spiral example
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
T = T(T.face &T.eye_DLC,:);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));
BW = logical(projectedAtlas1);                                             % atlas brain boundry binary mask 
[row,col] = find(BW);
brain_index = [col,row];  
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% Parameters
fs = 35;                    % Sampling frequency
low_freq_band = [0.05, 0.5]; % Phase-providing frequency (slow oscillation)
high_freq_band = [2, 8];    % Amplitude-modulated frequency (fast oscillation)
%% session info
kk = 3;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
% read svd components (U,V,t) from processed data folder
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
serverRoot = expPath(mn, td, en);
session_root = fullfile(data_folder,'spirals','svd',subfolder);
%%
[U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
% registration
load(fullfile(data_folder,'spirals','rf_tform',...
    [fname '_tform.mat']));  
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
serverRoot = expPath(mn, td, en);
v = VideoReader(fullfile(serverRoot,'face.mp4'));
frameID = 3030*70;
frame = read(v,frameID);
frame = rgb2gray(frame);
%%
v2 = VideoReader(fullfile(serverRoot,'eye.mp4'));
frame_eye = read(v2,frameID);
frame_eye = rgb2gray(frame_eye);
%%
hs= figure('Renderer', 'painters', 'Position', [100 100 600 300]);
ax1 = subplot(1,2,1);
frame1 = imrotate(frame,-90);
imagesc(frame1);
colormap(ax1,gray);
axis image;axis off; 
caxis([min(frame1(:)),max(frame1(:))*0.8]);

ax2 = subplot(1,2,2);
imagesc(frame_eye);
colormap(ax2,gray);
axis image;axis off; 
%%
print(hs, ['face_eye_frame.pdf'],'-dpdf', '-bestfit', '-painters');
%% pupil
tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
pupil_mean = readNPY(tFile1); 
pupil_mean2 = pupil_mean(1:2:end);
if numel(pupil_mean2)<size(t,1)
    pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
elseif numel(pupil_mean2)>size(t,1)
    pupil_mean2 = pupil_mean2(1:numel(t));
end 
%% rotaryEncoder
sigName = 'rotaryEncoder';
tlFile = fullfile(serverRoot, [sigName '.raw.npy']); 
pd = readNPY(tlFile);
fs1 = 35;
wh = wheel.correctCounterDiscont(pd);
tlFile = fullfile(serverRoot, [sigName '.timestamps_Timeline.npy']);
tlTimes = readNPY(tlFile);
tt = tsToT(tlTimes, numel(pd));  
fs1 = 1/mean(diff(tt));
wh = wheel.correctCounterDiscont(pd);
vel = wheel.computeVelocity2(wh, 0.01, fs1); 
vel2 = interp1(tt,vel,t);
%% motion energy
load(fullfile(data_folder,'spirals','spirals_index',...
[fname '_motion_energy.mat']));
image_energy2(isnan(image_energy2)) = 0;
if numel(image_energy2)<size(t,1)
    image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
elseif numel(image_energy2)>size(t,1)
    image_energy2 = image_energy2(1:numel(t));
end 
%%  filter ficial motion
% Design filters for frequency bands
% Low frequency (phase) - 4th order Butterworth
[b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
signal_low = filtfilt(b_low, a_low, image_energy2);    
% Extract phase from low frequency signal using Hilbert transform
phase_low = angle(hilbert(signal_low));

[b1, a1] = butter(4, 8/(fs/2));
image_energy3 = filtfilt(b1, a1, image_energy2);  
%% licks
tFile1 = fullfile(serverRoot, 'lick.raw.npy'); 
lick = readNPY(tFile1); 
lick2 = lick(1:2:end);
if numel(lick2)<size(t,1)
    lick2(numel(lick2)+1:size(t,1)) = 0;
elseif numel(lick2)>size(t,1)
    lick2 = lick2(1:numel(t));
end 
%%
fs = 35;
[b2, a2] = butter(4, 1/(fs/2),'high');
lick3 = filtfilt(b2, a2, lick2);  
[pks,locs] = findpeaks(lick3);
%% reward valve
tFile1 = fullfile(serverRoot, 'rewardValve.raw.npy'); 
reward = readNPY(tFile1); 
tlFile = fullfile(serverRoot, 'rewardValve.timestamps_Timeline.npy');
tlTimes = readNPY(tlFile);
tt = tsToT(tlTimes, numel(reward)); 
reward2 = interp1(tt,reward,t);
%%
figure;
plot(t,reward2)
[allPD, flipsUp, flipsDown] = schmittTimes(t,reward2, [0.5 1]);
%%
clear trace signal signal_high amplitude_high
[b_high, a_high] = butter(4, high_freq_band/(fs/2), 'bandpass');
for j = 1:7
    trace(:,j) = squeeze(Utransformed(pixel(j,1),pixel(j,2),1:50))'*V(1:50,:);
    signal(:,j) = double(trace(:,j))'./mimgtransformed(pixel(j,1),pixel(j,2));
    signal_high(:,j) = filtfilt(b_high, a_high, signal(:,j));
    amplitude_high(:,j) = abs(hilbert(signal_high(:,j)));
end
%% 
clear spiralsT filteredSpirals lia locb
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
%%
spiral_count = zeros(numel(t),1);
for jj = 1:numel(t)
    spiral_count(jj,1) = sum(filteredSpirals(:,5)>=jj-17 &filteredSpirals(:,5)<jj+17);
end
%%
h1= figure('Renderer', 'painters', 'Position', [100 100 600 400]);
subplot(1,1,1);
locs1 = locs(pks>3);
t1 = t(locs1);
t2 = t1(t1>3000 &t1<3060);
amp_mean = mean(amplitude_high(:,3:7),2);
amp_mean2 = smoothdata(amp_mean,"gaussian",35);

plot(t,spiral_count/10,'c');
% ylim([1.5,2]);
yticks([1.5:0.25:2]);

hold on;
plot(t,reward2/4,'g');
hold on;
plot(t,image_energy3/2000000+0.4,'r');
hold on;
plot(t,pupil_mean2/10-1.2,'color',[0.6,0.6,0.6]);
hold on;
scatter(t2,ones(numel(t2))*0.2,'|');
hold on;
% plot(t,lick3/80,'color',[0.6,0.6,0.6]);


for k = 3:7
    plot(t,signal(:,k)*30-1*k,'k');
    hold on;
end
plot([3000,3000],[-2,0],'r')

hold on;
plot(t,amp_mean2*100,'m');
hold on;
plot([3000,3000],[0,0.5],'r')

xlim([3020,3050]);
%%
h1= figure('Renderer', 'painters', 'Position', [100 100 600 400]);
locs1 = locs(pks>3);
t1 = t(locs1);
t2 = t1(t1>3000 &t1<3060);
amp_mean = mean(amplitude_high(:,3:7),2);
amp_mean2 = smoothdata(amp_mean,"gaussian",35);
subplot(4,1,1);
plot(t,spiral_count,'c');
% ylim([0,0.8]);
yticks([0:12:24]);
xlim([3020,3050]);

subplot(4,1,2);
plot(t,image_energy3/2000000+0.4,'r');
hold on;
plot(t,pupil_mean2/10-1.2,'color',[0.6,0.6,0.6]);
% plot(t,pupil_mean2/2-9,'color',[0.6,0.6,0.6]);
hold on;
scatter(t2,ones(numel(t2))*0.2,'|');
% hold on;
% plot(t,lick3/80,'color',[0.6,0.6,0.6]);
xlim([3020,3050]);

subplot(4,1,3);
for k = 3:7
    plot(t,signal(:,k)*100-2*k,'k');
    hold on;
end
plot([3020,3020],[-2,0],'r')
xlim([3020,3050]);

subplot(4,1,4);
plot(t,amp_mean2*100,'m');
hold on;
plot([3020,3020],[0,0.5],'r')
xlim([3020,3050]);
%%
print(h1, [mn '_example2.pdf'],'-dpdf', '-bestfit', '-painters');