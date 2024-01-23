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
spiral_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\fftn_permutation';
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
spirals_all = [];
frame_all = 0;
for kk = 1:session_total
    clear spiralsT
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%    
    dfolder = fullfile(folder,mn,td,num2str(en));
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(spiral_folder,[fname '.mat']));

    tformName = dir(fullfile(dfolder,'*tform.mat')).name;
    load(fullfile(dfolder,tformName));

    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,pwAll1(:,1),pwAll1(:,2));
    pwAll1(:,1:2) = round(spiralsT); 
    spirals_all = [spirals_all;pwAll1];
    frame_all = frame_all+frame_count;
end
%%
spirals_all(:,1:2) = round(spirals_all(:,1:2));
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
[lia,locb] = ismember(spirals_all(:,1:2),brain_index,'rows');
spirals_all = spirals_all(lia,:);
%%
hist_bin = 40;
[unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spirals_all,hist_bin);
save(['histogram_' num2str(hist_bin) 'pixels.mat'], 'unique_spirals');
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
spiral_count = size(spirals_all,1);
hist_bin = 40;
h = figure;
clear unique_spiral_unit
load(['histogram_' num2str(hist_bin) 'pixels.mat']); %load histogtam counts
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
ax1 = subplot(1,2,1);
scatter(unique_spirals(:,1),unique_spirals(:,2),3,unique_spirals_unit,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
colormap(ax1,hot);
% axis off; 
axis image;
xlim(ax1,[0,1140]);
ylim(ax1,[0,1320]);
xtick = 0:200:1000;
xtick_string = string(xtick);
xticks(xtick)
xticklabels(xtick_string)
ytick = 0:200:1200;
ytick_string = string(ytick);
yticks(ytick)
yticklabels(ytick_string)
% caxis([0,130]);
cb1 = colorbar;
cb_size1 =  get(cb1,'Position'); %gets the positon and size of the color bar
set(cb1,'Position',[cb_size1(1)+0.07  cb_size1(2) cb_size1(3) cb_size1(4)/2]);% To change size
% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
hold on;
scatter(points(:,1),points(:,2),8,'g','filled');
% interp histgram counts 
F = scatteredInterpolant(unique_spirals(:,1),unique_spirals(:,2),unique_spirals_unit );
xx = 1:1140; yy = 1:1320;
[xxq,yyq] = meshgrid(xx,yy);
vq1 = F(xxq,yyq);
for i = 1:size(points,1)
    count_sample(i,1) = vq1(points(i,2),points(i,1));
end
subplot(1,2,2);
count_sample(count_sample<0) = 0;
% plot(points(:,1),count_sample);
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
plot(points_line,count_sample);
xlim([0,800]);
% ylim([0,4.5]);
xticks([0 100 200 300 400 500 600 700 800]);
xticklabels({'0','1','2','3','4','5','6','7','8'});
xlabel('Distance (mm)');
ylabel('Sprials/(mm^2*s)');
print(h, 'spiral_fftn_histogram_with_40binsize', '-dpdf', '-bestfit', '-painters');
%%
hist_bin = 40;
for radius = 10:10:100
    clear spirals_temp
    spirals_temp = spirals_all(spirals_all(:,3)==radius,:);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spirals_temp,hist_bin);
    save(['histogram_fftn_radius' num2str(radius) '.mat'], 'unique_spirals');
end
%%
