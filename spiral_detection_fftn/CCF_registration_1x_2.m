githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
% addpath(genpath('C:\Users\Steinmetz lab\Documents\git\AP_scripts_cortexlab-master'))
addpath(genpath(fullfile(githubDir, 'wheelAnalysis')))
addpath(genpath(fullfile(githubDir, 'PatternDetection')))
addpath(genpath(fullfile(githubDir, 'spiralDetection')));
addpath(genpath(fullfile(githubDir, 'allenCCF')));
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
BW = logical(projectedAtlas1);
%% get original and predicted dV by kernel regression.
% codefolder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\PowerMap';
% T = readtable(fullfile(codefolder, 'waterRewardMice.xlsx'));
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
%%
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
%%
% for kk = [1,4,5,8,9]
for kk = 11:17
    %%
    clearvars -except projectedAtlas1 projectedTemplate1 T kk indx gfolder scale coords index_all rf_folder BW
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    serverRoot = expPath(mn, td, en);
    %%
    fname = [mn '_' tdb '_' num2str(en) '_tform.mat'];
    dfolder = fullfile(gfolder,mn,td,num2str(en));
    ffname =fullfile(dfolder,fname);
    %%
    if ~(exist(ffname, 'file') == 2)
        mimg = readNPY(fullfile(serverRoot, 'blue', ['meanImage.npy']));
        mimg1 = mimg/max(mimg(:));
        mimg1 = mimg1(1:scale:end,1:scale:end);
        clear cpstruct
        projectedTemplate1 = projectedTemplate1/max(projectedTemplate1(:));
        %%
        h = cpselect(mimg1,projectedTemplate1);
        while ~exist('cpstruct','var')
            pause(2)
        end
        [movingPoints,fixedPoints] = cpstruct2pairs(cpstruct); 
        tform = fitgeotrans(movingPoints,fixedPoints,'affine');    
        save(ffname,'tform');
        %%
        sizeTemplate = size(projectedTemplate1);
        mimgtransformed = imwarp(mimg1,tform,'OutputView',imref2d(sizeTemplate));
        figure
        ax1 = subplot(1,2,1); 
        im = imagesc(ax1, mimgtransformed);
        set(im, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
        colormap(ax1,'gray');
        hold on; 
        scale1 = 1;
        overlayOutlines(coords,scale1,'k');
        % set(gca,'Ydir','reverse');
        axis off; axis image;
        %%      
        rfname = [mn '_' tdb '_' num2str(en) '_signMap'];
        % load(fullfile(rf_folder,dir(fullfile(rf_folder, [mn '*'])).name));
        load(fullfile(dfolder,[rfname, '.mat']));
        %%
        ax2 = subplot(1,2,2); 
        mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
        signMaptransformed = imwarp(signMap,tform,'OutputView',imref2d(size(projectedTemplate1)));   
        im2 = imagesc(ax2,signMaptransformed);
        % set(ax2, 'Color', 'w');
        set(ax2, 'XTickLabel', '', 'YTickLabel', '');
        axis image; axis off;
        box off
        % colormap(ax0,'gray')
        set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
        colormap(ax2,colormap_RedWhiteBlue)
        hold on;
        overlayOutlines(coords,scale1,'k');
    end
end