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
% only select cortex in the atlas
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'ara_tools-master')))
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
areaPath{1} = '/997/8/567/688/695/315/453/322/'; % SSp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
scale = 1;
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_right = BW_empty; 
BW_right(indexright) =1;
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
[row,col] = find(BW_right);
brain_index = [col,row];
%%
control_folder  = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\control';
permute_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\sessions\permute';
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
count_sample_control = zeros(15,10);
count_sample_permute = zeros(15,10);
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
    count = 1;
    for radius = 10:10:100
        clear unique_spirals_unit_control unique_spirals_unit_permute F_control F_permute vq1_control vq1_permute
        spirals_control = control.spiral_density{count};
        spirals_permute = permute.spiral_density{count};
        [lia,locb] = ismember(spirals_control(:,1:2),brain_index,'rows');
        spirals_control = spirals_control(lia,:);
        [lia1,locb1] = ismember(spirals_permute(:,1:2),brain_index,'rows');
        spirals_permute = spirals_permute(lia1,:);
        unique_spirals_unit_control = spirals_control(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = spirals_permute(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        % interp histgram counts 
        count_sample_control(kk,count) = max(unique_spirals_unit_control(:));
        count_sample_permute(kk,count) = max(unique_spirals_unit_permute(:));
        count = count+1;
    end
end
%%
h = figure;
for radius = 1:10
    subplot(1,10,radius)
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    scatter(1,control_all,4,'k','filled');
    hold on;
    scatter(2,permute_all,4,'k','filled');
    hold on;
    x_control = ones(size(control_all)); x_permute = ones(size(control_all))*2;
    plot([x_control,x_permute]',[control_all,permute_all]','k');
    ylim([0,2.5]);
end
print(h, 'density_significance_across_session', '-dpdf', '-bestfit', '-painters');
%%
for radius = 1:10
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    [ha(radius),pa(radius)] = ttest(control_all,permute_all);
end
%% percentage change
for radius = 1:10
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    percentage_change(:,radius) = (control_all-permute_all)./permute_all;
end
%%
for i = 1:10
    percentage_temp = percentage_change(:,i);
    [h3(i),p3(i)] = ttest(percentage_temp);
end
%%
mean_percentage_change = mean(percentage_change,1);
std_percentage_change = std(percentage_change,[],1);
sem_percentage_change = std_percentage_change/sqrt(15);
%%
x_point = normrnd(0,0.1,[1,15]);
%%
scatter_x = repmat([1:10],15,1);
scatter_x = scatter_x+repmat(x_point',1,10);
%%
h = figure;
scatter(scatter_x,percentage_change*100,6,'k','filled');
hold on;
x = 1:10;
bar(x,mean_percentage_change*100,'FaceColor','None');
hold on
er = errorbar(x,mean_percentage_change*100,sem_percentage_change*100);    
er.Color = [0,0,0];                         
er.LineStyle = 'none';  
hold off
ylim([-100 500]);
xlabel('Radius');
ylabel('Percentage change (%)');
% print(h, 'density_percentage_change_across_session', '-dpdf', '-bestfit', '-painters');