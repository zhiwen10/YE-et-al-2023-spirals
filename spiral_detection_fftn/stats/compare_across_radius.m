% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
BW = logical(projectedAtlas1);
%%
session_all = find(T.use);
session_total = numel(session_all);
%%
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
control_folder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\control';
permute_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\permute';
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
for kk = 1:session_total
    clear spiralsT
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    control = load(fullfile(control_folder,[mn '_density.mat']));
    permute = load(fullfile(permute_folder,[mn '_density.mat']));
    h = figure('Renderer', 'painters', 'Position', [100 100 900 900]);
    count = 1;
    for radius = 10:10:50
        %%
        unique_spirals_unit_control = control.spiral_density{count}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = permute.spiral_density{count}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        %%
        cmax = max([unique_spirals_unit_control;unique_spirals_unit_permute]);
        %%
        ax1 = subplot(5,4,4*count-3);
        scatter(control.spiral_density{count}(:,1),control.spiral_density{count}(:,2),3,unique_spirals_unit_control,'filled');
        hold on;
        scale2=1;
        overlayOutlines(coords,scale2);
        set(gca,'Ydir','reverse')
        colormap(ax1,hot);
        axis image; axis off;
        xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
        caxis([0,cmax]);
        xtick = 0:200:1000; xtick_string = string(xtick);
        xticks(xtick); xticklabels(xtick_string);
        ytick = 0:200:1200; ytick_string = string(ytick);
        yticks(ytick); yticklabels(ytick_string);
        colorbar;    
        ax2 = subplot(5,4,4*count-2);
        scatter(permute.spiral_density{count}(:,1),permute.spiral_density{count}(:,2),3,unique_spirals_unit_permute,'filled');
        hold on;
        scale2=1;
        overlayOutlines(coords,scale2);
        set(gca,'Ydir','reverse')
        colormap(ax2,hot);
        axis image; axis off;
        xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
        caxis([0,cmax])
        colorbar;
        count = count+1;
    end
    %%
    count = 1;
    for radius = 60:10:100
        %%
        unique_spirals_unit_control = control.spiral_density{count+5}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = permute.spiral_density{count+5}(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        %%
        cmax = max([unique_spirals_unit_control;unique_spirals_unit_permute]);
        %%
        ax1 = subplot(5,4,4*count-1);
        scatter(control.spiral_density{count+5}(:,1),control.spiral_density{count+5}(:,2),3,unique_spirals_unit_control,'filled');
        hold on;
        scale2=1;
        overlayOutlines(coords,scale2);
        set(gca,'Ydir','reverse')
        colormap(ax1,hot);
        axis image; axis off;
        xlim(ax1,[0,1140]); ylim(ax1,[0,1320]);
        caxis([0,cmax]);
        xtick = 0:200:1000; xtick_string = string(xtick);
        xticks(xtick); xticklabels(xtick_string);
        ytick = 0:200:1200; ytick_string = string(ytick);
        yticks(ytick); yticklabels(ytick_string);
        colorbar;    
        ax2 = subplot(5,4,4*count);
        scatter(permute.spiral_density{count+5}(:,1),permute.spiral_density{count+5}(:,2),3,unique_spirals_unit_permute,'filled');
        hold on;
        scale2=1;
        overlayOutlines(coords,scale2);
        set(gca,'Ydir','reverse')
        colormap(ax2,hot);
        axis image; axis off;
        xlim(ax2,[0,1140]); ylim(ax2,[0,1320]);
        caxis([0,cmax])
        colorbar;
        count = count+1;
    end
    print(h, [mn '_radius'], '-dpdf', '-bestfit', '-painters');
end
