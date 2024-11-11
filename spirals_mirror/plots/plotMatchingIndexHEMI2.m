function h3i = plotMatchingIndexHEMI2(data_folder,save_folder)
%%
load(fullfile(data_folder,'spirals_mirror\matching_index',...
    'hemi_weights_8points_allsessions.mat'));
load(fullfile(data_folder,'spirals_mirror\matching_index',...
    'axon_intensity_all_hemi.mat'));
gcamp_mean = mean(TheColorImage_all,4);
%%
for i = 1:8
    intensity_all1(:,:,i) = imresize(squeeze(intensity_all(:,:,i)),...
        size(TheColorImage_all,[1,2]));
end
%%
dot_real = 0;
for kkk = 1:8
    gcamp = squeeze(gcamp_mean(:,:,kkk));
    axon = squeeze(intensity_all1(:,:,kkk));
    dot_temp = dot(gcamp(:),axon(:));
    dot_real = dot_real+dot_temp;
end
dot_real_all2 = dot_real;
%
for i = 1:1000
    interation = randperm(8);
    intensity_all2 = intensity_all1(:,:,interation);
    dot_perm = 0;
    for kkk = 1:8
        gcamp = squeeze(gcamp_mean(:,:,kkk));
        axon = squeeze(intensity_all2(:,:,kkk));
        dot_temp = dot(gcamp(:),axon(:));
        dot_perm = dot_perm+dot_temp;
    end
    dot_perm_all2(i,1) = dot_perm;
end
[h2,p2] = ttest2(dot_real_all2, dot_perm_all2);
%%
dot_max = 0;
for kkk = 1:8
    axon = squeeze(intensity_all1(:,:,kkk));
    dot_temp = dot(axon(:),axon(:));
    dot_max = dot_max+dot_temp;
end
%% normalize between 0 and 1
dot_all2 = [dot_real_all2; dot_perm_all2];
% dot_all2 = dot_all2-min(dot_all2(:));
% dot_all2 = dot_all2/max(dot_all2(:));
dot_all2 = dot_all2./dot_max;
%% significance test
mean_indx = mean(dot_all2(2:1001));
std_indx = std(dot_all2(2:1001));
sem_indx  =std_indx./sqrt(1000);
[h2,p2] = ttest2(dot_all2(1), dot_all2(2:1001));
%%
h3i = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
% edges = [0:0.05:1];
% histogram(dot_all2(2:end),edges)
histogram(dot_all2(2:end))
hold on;
xline(dot_all2(1),'r')
% xlim([0,1.2]);
%%
print(h3i, fullfile(save_folder,'Fig3i_matching_index_hemi.pdf'),...
    '-dpdf', '-bestfit', '-painters');