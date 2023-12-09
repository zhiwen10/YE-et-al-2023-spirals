githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')));
addpath(genpath(fullfile(githubdir2, 'Pipelines')));
addpath(genpath(fullfile(githubdir2, 'widefield')));
addpath(genpath(fullfile(githubdir2, 'npy-matlab')));
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection');
% addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox'))) %https://github.com/BrainDynamicsUSYD/NeuroPattToolbox
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
T = readtable('spiralSessions.xlsx');
session_all = find(T.use);
session_total = numel(session_all);
T1 = T(session_all,:);
id = [1,2,3,7,10];
%%
for kk = numel(id)
    %%
    mn = T1.MouseID{id(kk)};
    tda = T1.date(id(kk));
    en = T1.folder(id(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% SVD, plot example trace, overlay with pupil and stim time
    serverRoot = expPath(mn, td, en);
    mimg = readNPY(fullfile(serverRoot, 'blue','meanImage.npy'));
    params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
    params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
    params.padding = 2*params.halfpadding;
    %% apply mask, this helps speed up spiral detection later
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    %%
    figure; 
    ax1 = imagesc(mimg2);
    roi = drawpolygon;
    fname1 = [mn '_' tdb '_' num2str(en) '_roi'];
    save(fname1,'roi');
end