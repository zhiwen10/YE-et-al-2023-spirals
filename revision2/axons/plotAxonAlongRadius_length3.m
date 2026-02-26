function [h2ac] = plotAxonAlongRadius_length3(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells2.csv'));
%% examples
h2ac = figure('Renderer', 'painters', 'Position', [100 100 500 200]);
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
center = [244,542]; % SSp-un
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
soma_ref = soma_center1-center;
distance = vecnorm(soma_ref,2,2)*0.01;
axon_length = T1.axon_length*0.01;
cellN = sum(distance<=4);
distance = distance(distance<=4);
axon_length = axon_length(distance<=4);

ax4 = subplot(1,2,1);
scatter(distance,axon_length,6,'k');
hold on;
radius_bins = 0:0.25:4.5;
pc_ratio_bins = {};
for i = 1:length(radius_bins)-1
    pc_ratio_index = (distance>radius_bins(i) & distance<radius_bins(i+1));
    pc_ratio_temp = axon_length(pc_ratio_index);
    pc_ratio_bins{i} = pc_ratio_temp;
end
pc_ratio_bins_mean = cellfun(@mean, pc_ratio_bins);
pc_ratio_bins_std = cellfun(@std, pc_ratio_bins);
pc_ratio_bins_count = cellfun(@length,pc_ratio_bins);
pc_ratio_bins_sem = pc_ratio_bins_std./sqrt(pc_ratio_bins_count);
radius_bins2 = radius_bins(2:end);
index = (pc_ratio_bins_count>1);
errorbar(radius_bins2(index),pc_ratio_bins_mean(index), pc_ratio_bins_sem(index),'lineWidth',2);
ylim([0,5]);
box off;
xlabel('Soma distance to SSp center (mm)');
ylabel('Axon length along PC1 (mm)');

ax5 = subplot(1,2,2);
clear coef2
[coef1,p1] = corrcoef(distance, axon_length);
coef1 = coef1(1,2);
data_n = length(distance);
for i = 1:1000
    pa = randperm(data_n);
    distance1 = distance(pa);
    corr_temp = corrcoef(distance1,axon_length);
    coef2(i,1) = corr_temp(1,2);
end
histogram(coef2, 'BinEdges',[-0.2:0.02:0.2]);
hold on;
xline(coef1,'r--');
xlim([-0.2,0.6]);
box off;
xlabel('Correlation Coefficient');
ylabel('Permutation Counts');

[h,p] = ttest2(coef2,coef1);
%%
print(h2ac, fullfile(save_folder,'Fig2ac_axon_bias_all2_axon_length4.pdf'),...
    '-dpdf', '-bestfit', '-painters');