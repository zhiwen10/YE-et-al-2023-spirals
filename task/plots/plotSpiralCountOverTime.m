function hs15k = plotSpiralCountOverTime(data_folder,save_folder)
%% load atlas brain horizontal projection and outline 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
pix_sum = sum(BW(:));
area_cortex = pix_sum*pixArea;
%%
load(fullfile(data_folder,'task','spirals','task_spiral_count_by_radius_time'));
spiral_count_sum_all = spiral_count_sum_all./area_cortex;
color1 = {'k','g','r'};
hs15k = figure('Renderer', 'painters', 'Position', [100 100 900 600]);
for i = 1:3
    spiral_count_sum_all1 = squeeze(spiral_count_sum_all(i,:,:,:));
    spiral_count_sum_mean = squeeze(mean(spiral_count_sum_all1,1));
    spiral_count_sum_sem = squeeze(std(spiral_count_sum_all1,[],1))./sqrt(4);
    for j = 1:7
        subplot(3,7,j+(i-1)*7);
        spiral_count_sum = spiral_count_sum_mean(j,:);
        spiral_count_sem = spiral_count_sum_sem(j,:);
        shadedErrorBar(1:141, spiral_count_sum, spiral_count_sem, 'lineprops', color1{i});
        xline(71,'--k');
        xticks([1,35,70,105,141]);
        xticklabels({'-2','-1','0','1','2'});
        ylim([0,0.08]);
    end
end
%%
print(hs15k, fullfile(save_folder,'FigS15k_task_spiral_count_over_time'), ...
    '-dpdf', '-bestfit', '-painters');