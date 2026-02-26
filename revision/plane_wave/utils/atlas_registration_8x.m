githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\AP_scripts_cortexlab-master'))
addpath(genpath(fullfile(githubDir, 'wheelAnalysis')))
addpath(genpath(fullfile(githubDir, 'PatternDetection')))
%%
save_folder = 'E:\spiral_data_share\data\spirals\rf_tform_8x';

%% load 10um horizontal atlas and outline
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\projectedOutlineAtlas.mat')
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'ara_tools-master')))
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
% tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
% av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%% only select cortex in the atlas
projectedAtlas1 = projectedAtlas;
projectedTemplate1 = projectedTemplate;
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/512/'); %CB
spath3 = startsWith(spath,'/997/8/343/'); %BS
idFilt2 = st.index(spath2);
idFilt3 = st.index(spath3);
idFilt = [idFilt2;idFilt3];
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(Lia) = 0; 
projectedTemplate1(Lia) = 0;
projectedAtlas2 = projectedAtlas;
spath4 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath4);
st1  = st(spath4,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas2(~Lia) = 0; 
projectedTemplate2(~Lia) = 0;
projectedAtlas1(700:end,:) = projectedAtlas2(700:end,:);
%%
scale = 8;
projectedTemplate1 = projectedTemplate1(1:scale:end,1:scale:end);
%%
for kk = 1:15
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals\svd',subfolder);
    mimg = readNPY(fullfile(session_root, 'meanImage.npy'));
    %% outline binary image
    mimg1 = mimg/max(mimg(:));
    mimg2 = mimg1(1:scale:end,1:scale:end);
    clear cpstruct
    projectedTemplate1 = projectedTemplate1/max(projectedTemplate1(:));
    h = cpselect(mimg2,projectedTemplate1);

    while ~exist('cpstruct','var')
        pause(2)
    end
    [movingPoints,fixedPoints] = cpstruct2pairs(cpstruct); 
    tform = fitgeotrans(movingPoints,fixedPoints,'affine');
    %% save tform 
    fname1 = [fname '_tform_8x'];
    save(fullfile(save_folder,fname1),'tform');
    %%
    sizeTemplate = size(projectedTemplate1);
    mimgtransformed = imwarp(mimg2,tform,'OutputView',imref2d(sizeTemplate));
    %%
    figure
    ax2 = subplot(1,1,1); 
    % im = imagesc(signMaptransformed);
    % set(im, 'AlphaData', mimgtransformed, 'AlphaDataMapping', 'scaled');
    im = imagesc(mimgtransformed);
    set(im, 'AlphaData', logical(mimgtransformed), 'AlphaDataMapping', 'scaled');
    % set(ax2, 'Color', 'k');
    % set(ax2, 'XTickLabel', '', 'YTickLabel', '');
    axis image
    box off; axis off;
    % colormap(colormap_RedWhiteBlue)
    hold on; 
    for q = 1:numel(coords) % coords is from ctxOutlines.mat 
        % these are in 10um voxel coordinates, so first convert to mm, then to
        % pixels 
        cx = coords(q).x/8;
        cy = coords(q).y/8;    
        coordsX(q).x = cx;
        coordsX(q).y = cy;
        plot(cx,cy, 'LineWidth', 1.0, 'Color', 'r');
        hold on;
    end
    %%
    close all force
end