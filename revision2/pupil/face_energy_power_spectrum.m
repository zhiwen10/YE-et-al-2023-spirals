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
T = T(logical(T.face),:);
% T = T([1,2,6:13],:);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
freq_value = [0.05,0.1,0.2,0.5,1,2,4,6,8,10];
log_freq_value = log10(freq_value);
freqN = 350;
%% Parameters
fs = 35;                    % Sampling frequency
figure;
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    serverRoot = expPath(mn, td, en);
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    %%
    load(fullfile(data_folder,'spirals','spirals_index',...
    [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    signal = image_energy2;
    %%
%     tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
%     pupil_mean = readNPY(tFile1); 
%     pupil_mean2 = pupil_mean(1:2:end);
%     signal = pupil_mean2;
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
%     signal = vel2;
    %%
    epochLength = 20*35;
    epochs = floor(numel(signal)./epochLength);
    for epoch_i = 1:epochs
        firstFrame(epoch_i,1) = 1+(epoch_i-1)*epochLength;
        lastFrame(epoch_i,1) = epoch_i*epochLength;
        traceEpoch(epoch_i,:) = signal(firstFrame(epoch_i):lastFrame(epoch_i));
    end
    [freq1,psdx,~] = fft_spectrum(traceEpoch);
    psdx_mean = mean(psdx,2);
    %%
    plot(log10(freq1(2:freqN)),log10(psdx_mean(2:freqN)),'color','k');
    hold on;
    xlim([log_freq_value(1),log_freq_value(end)]);
    xticks(log_freq_value(1:end));
    xticklabels({'0.05','0.1','0.2','0.5','1','2','4','6','8','10'});
    xlabel('log10(Frequency)');
    ylabel('log10(Power) (df/f^2)');
end
