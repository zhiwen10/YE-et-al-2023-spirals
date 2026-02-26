%% pupil_amp vs rate correlation coefficient
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
%% Parameters
fs = 35;                                                                   % Sampling frequency
%%
for kk = 1:size(T,1)
    %%
    clear spiralsT spiral_count amp_mean2 amp_mean;
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
    Utransformed = imwarp(U,tform,'OutputView',imref2d([1320,1140]));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d([1320,1140]));
    %%
    tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
    pupil_mean = readNPY(tFile1); 
    pupil_mean2 = pupil_mean(1:2:end);
    if numel(pupil_mean2)<size(t,1)
        pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
    elseif numel(pupil_mean2)>size(t,1)
        pupil_mean2 = pupil_mean2(1:numel(t));
    end 
    pupil = pupil_mean2;
    %% Design filters for frequency bands
    % Low frequency (phase) - 4th order Butterworth
    [b_low, a_low] = butter(4, low_freq_band/(fs/2), 'bandpass');
    pupil_low = filtfilt(b_low, a_low, pupil);    
    % pupil_low_amp = abs(hilbert(pupil_low));
    %%
    [b_high, a_high] = butter(4, high_freq_band/(fs/2), 'bandpass');
    amplitude_high = [];
    for j = 1:7
        trace = squeeze(Utransformed(pixel(j,1),pixel(j,2),1:50))'*V(1:50,:);
        signal = double(trace)'./mimgtransformed(pixel(j,1),pixel(j,2));
        signal_high = filtfilt(b_high, a_high, signal);
        amplitude_high(:,j) = abs(hilbert(signal_high));
    end
    amp_mean = mean(amplitude_high(:,3:7),2);
    % amp_mean2 = smoothdata(amp_mean,"gaussian",35);
    %% only use spiral sequences with at least 2 consecutive frames
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);                                          % sprial length > = 2 frames
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    spiral_count = zeros(numel(t),1);
    for jj = 1:numel(t)
        spiral_count(jj,1) = sum(filteredSpirals(:,5)>=jj-17 &filteredSpirals(:,5)<jj+17);
    end
    %%
    R1 = corrcoef(pupil_low,spiral_count);
    R2 = corrcoef(amp_mean,spiral_count);
    r1(kk,1) = R1(1,2);
    r2(kk,1) = R2(1,2);
end
%%
[h,p]= ttest(r1,r2);
%%
h1 = figure('Position', [100, 100, 300, 400]);
scatter(ones(length(r1)),r1,8,'k');
hold on;
scatter(ones(length(r2))*2,r2,8,'r');
plot([ones(length(r1),1)*1,ones(length(r2),1)*2]',...
    [r1,r2]','color','k');
xticks([1,2]);
xticklabels({'Pupil & rates','2-8Hz amp & rates'});
ylabel('Correlation coefficent');
text(1.5,0,['p = ' num2str(p)]);
print(h1, 'coeff','-dpdf', '-bestfit', '-painters');