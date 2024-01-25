githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%% cortex surface outline
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
%%
prediction_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\prediction';
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\roi';
%%
scale = 1;
scale3 = 5/scale;
hemi = [];
h1 = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);
count1 = 1;
for areai = [1,2,4]
    current_T = T(contains(T.Area,area(1,areai)),:);
    for kk = 1:size(current_T)
        pos = kk+(count1-1)*12;
        ax2 = subplottight(3,12,pos);
        ops = get_session_info(current_T,kk);
        fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
        fname_predict = fullfile(prediction_folder,[fname '_dv_predict.mat']);
        load(fname_predict,'explained_var_all');
        load(fullfile(roi_folder,[fname '_roi.mat']));
        explained_var_all(explained_var_all<0) = 0;
        var_mean = squeeze(mean(explained_var_all,3));
        fname1 = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x'];
        regist_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration_4x';
        load(fullfile(regist_folder, fname1));
        mean_var_t = imwarp(var_mean,tform,'OutputView',imref2d(size(projectedTemplate1)));
        im2 = imagesc(mean_var_t);
        set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
        caxis([0,0.8]);
        axis image; axis off;
        % cmap1 = flipud(inferno);
        cmap1 = inferno;
        colormap(ax2,cmap1);
        % colorbar;
        hold on;
        plotOutline(maskPath(1:3),st,atlas1,hemi,scale3);
        plotOutline(maskPath(4),st,atlas1,[],scale3);
        plotOutline(maskPath(5),st,atlas1,[],scale3);
        plotOutline(maskPath(6:11),st,atlas1,[],scale3);
        text(0,50,ops.mn,'Interpreter', 'none');
    end
    count1 = count1+1;
end
ax3 = subplottight(3,12,32);
hold off;
im2 = imagesc(mean_var_t);
set(im2, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
caxis([0,0.8]);
axis image; axis off;
% cmap1 = flipud(inferno);
cmap1 = inferno;
colormap(ax3,cmap1);
hold on;
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
cb = colorbar;
cb.Ticks = [0,0.4,0.8];
cb.TickLabels = string([0,0.4,0.8]);
print(h1, 'all_variance_map-inverted', '-dpdf', '-bestfit', '-painters');