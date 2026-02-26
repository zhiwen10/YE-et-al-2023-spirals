function [h5e] = plotSpiralRatePeriStimulus(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
scale = 1;
[row,col] = find(BW);
brain_index = [col,row];
%%
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
load(fullfile(data_folder,'whisker','spirals_peri_stim','Whisker_spirals_peri_stimulus.mat'));
%%
ylim2 = 12;
mouseN = 5;
t1 = -1:1/35:1;
h5e = figure('Renderer', 'painters', 'Position', [100 100 400 250]);
subplot(1,2,1);
spiral_count_sum_left1 = squeeze(sum(spiral_count_sum_left,2));
spiral_mean_left = mean(spiral_count_sum_left1,1);
spiral_sem_left = std(spiral_count_sum_left1,1)./sqrt(mouseN);
shadedErrorBar(t1, spiral_mean_left, spiral_sem_left, 'lineprops', 'g');
xline(t1(36),'k--');
xline(t1(39),'k--');
xline(t1(46),'k--');
ylim([0,ylim2]);
ylabel('Rotating waves/s');

subplot(1,2,2);
spiral_count_sum_right1 = squeeze(sum(spiral_count_sum_right,2));
spiral_mean_right = mean(spiral_count_sum_right1,1);
spiral_sem_right = std(spiral_count_sum_right1,1)./sqrt(mouseN);
shadedErrorBar(t1, spiral_mean_right, spiral_sem_right, 'lineprops', 'g');
xline(t1(36),'k--');
xline(t1(39),'k--');
xline(t1(46),'k--');
ylim([0,ylim2]);
ylabel('Rotating waves/s');

%%
print(h5e, fullfile(save_folder,'Fig5e_WhiskerSpiralsPeriStimulus.pdf'), '-dpdf', '-bestfit', '-painters');
