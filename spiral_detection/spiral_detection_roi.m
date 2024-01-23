githubdir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir, 'spikes'))); % https://github.com/cortex-lab/spikes
addpath(genpath(fullfile(githubdir, 'Pipelines'))); % https://github.com/SteinmetzLab/Pipelines
addpath(genpath(fullfile(githubdir, 'widefield'))); % https://github.com/cortex-lab/widefield
addpath(genpath(fullfile(githubdir, 'npy-matlab'))); % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals'))); 
%%
T1 = readtable('spiralSessions3.xlsx');
id = 1:15;
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