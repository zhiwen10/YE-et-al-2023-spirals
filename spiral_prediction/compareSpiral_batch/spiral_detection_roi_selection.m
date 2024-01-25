githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Colormaps'));
%% get original and predicted dV by kernel regression.
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
T = readtable('session_list_sorted2.csv');
T1 = T([contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1],:);
%%
params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35; % frame sampling rate
params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
params.padding = 2*params.halfpadding;
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\roi';
roi_exist = zeros(size(T1,1),1);
count1 = 1;
for kk = 28
% for kk = 1:size(T1,1)
    clearvars -except T1 kk params roi_folder roi_exist count1
    ops = get_session_info(T1,kk);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    serverRoot = expPath(ops.mn, ops.td, ops.en);   
    mimg = readNPY(fullfile(serverRoot, 'blue','meanImage.npy'));
    %% apply mask, this helps speed up spiral detection later
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    %%
    fname1 = [fname '_roi.mat'];
    filename = fullfile(roi_folder,fname1);
    if exist(filename, 'file') == 2
        roi_exist(kk) = 1;
    else
        roi_exist(kk) = 0;
        figure; 
        ax1 = imagesc(mimg2);
        roi = drawpolygon;
        save(fname1,'roi');
        close all;
    end
    count1 = count1+1;
end

