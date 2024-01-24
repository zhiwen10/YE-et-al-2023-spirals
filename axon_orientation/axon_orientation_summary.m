%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
% addpath(genpath(fullfile(githubDir, 'FMAToolbox')))
addpath(genpath(fullfile(githubDir, 'YE-et-al-2023-spirals'))); 
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
%%
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
load('C:\Users\Steinmetz lab\Documents\git\BIL\all_cell_with_parents.mat');
%%
axon_region_label = {'VIS','RSP','AUD','TEa','VISa','VISC','VISrl','SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un','SSs'};
axon_st_cortex_indx = get_region_index_list(axon_region_label,st); % get all index values that contain "Isocortex" in st table
axon_region_label2 = {'CTX','cc'};
axon_st_cortex_indx2 = get_region_index_list(axon_region_label2,st); % get all index values that contain "Isocortex" in st table
%%
AUD_label = {'AUD','TEa','VISC'};
AUD_indx = get_region_index_list(AUD_label,st); 
SSs_label = {'SSs'};
SSs_indx = get_region_index_list(SSs_label,st); 
SSp_label = {'VISrl','SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un'};
SSp_indx = get_region_index_list(SSp_label,st); 
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%% get horizontal brain outline 
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
% cortex surface outline
root1 = {'/997/'};
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
%% get soma and axons by area
% axon_st_cortex_indx: sensory area only
% axon_st_cortex_indx2: entire cortex
regions_all1 = {'VIS','RSP','AUD','TEa','VISa','VISC','VISrl','SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un','SSs'};
axon_indx{3} = AUD_indx;
axon_indx{4} = AUD_indx;
for kk = [1,2,5,6,7]
    axon_indx{kk} = axon_st_cortex_indx;
end
for kk = 8:14
    axon_indx{kk} = SSp_indx;
end
axon_indx{15} = SSs_indx;
%%
hemi_terminal = 1;
hemi_axon_all = 2;
axon_points_all2 = {};
soma_all2 = {};
axon_terminal_all2 = {};
for k = 1:numel(regions_all1)
% for k = 2
    iregion = regions_all1(k);
    st_region_indx = get_region_index_list(iregion,st); % get all index values that contain "SSp-ll" in st table
    [cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
    % only get terminals in the left sensory cortex
    [~,axon_terminal_all1] = get_axon_branching_by_neuron(cell_id,allCoords,av,axon_indx{k}, hemi_terminal);
    % get axon points for the left hemispheres in the entire cortex 
    [soma_all1,axon_points_all1,dendrite_all1] = get_axon_points_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx2, hemi_terminal);
    axon_points_all2{k} = axon_points_all1;
    soma_all2{k} = soma_all1;
    dendrite_all2{k} = dendrite_all1;
    axon_terminal_all2{k} = axon_terminal_all1;
    labels{k} = repmat(string(iregion),numel(soma_all1),1);
end
%%
soma_all3 = soma_all2;
axon_terminal_all3 = axon_terminal_all2;
soma_all3 = cat(2, soma_all3{:});
axon_terminal_all3 = cat(2,axon_terminal_all3{:});
labels_all = cat(1,labels{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
labels_all = labels_all(a);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
T = get_axon_bias_table(soma_center,vector_all,angle1,polarity,labels_all);
%%
clear T1
T1 = T;
% T1.pc_ratio = T1.pc_ratio/200;
% h = plot_distribution(T1,coords);
h = plot_distribution_outline(T1,coords,st);
%%
print(h, 'all_axon_bias_version8_color', '-dpdf', '-bestfit', '-painters');
%%
figure;
ax4 = subplot(1,1,1);
[ax4,angle_diff_real,angle_diff_perm,edges,N1,N_rp] = plot_angle_histogram_and_stats(ax4,T1,coords);
%%
% a1 = deg2rad(90); a2 = a1*ones(100,1);
% [h mu ul ll] = circ_mtest(a2,deg2rad(90));
[h mu ul ll] = circ_mtest(angle_diff_real,deg2rad(90));
mean_rad = [mu,ll,ul];
mean_deg = rad2deg(mean_rad);
pval = circ_wwtest(angle_diff_real,angle_diff_perm); % bad
for i = 1:100
    [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end