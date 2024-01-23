%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
load('all_cell_with_parents.mat');
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
axon_points_all3 = axon_points_all2;
soma_all3 = cat(2, soma_all3{:});
axon_terminal_all3 = cat(2,axon_terminal_all3{:});
axon_points_all3 = cat(2,axon_points_all3{:});
labels_all = cat(1,labels{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
labels_all = labels_all(a);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
axon_points_all3 = axon_points_all3(a);
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
T = get_axon_bias_table(soma_center,vector_all,angle1,polarity,labels_all);
%%
clear a1 T1
a1 = not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"])));
T1 = T(a1,:);
h = plot_distribution(T1,coords);
%%
x1 = linspace(-pi,pi,180);
color2 = colorcet('C06','N',180);
color3 = interp1(x1,color2,T1.soma_angle);
center = [244,542]; % SSp-un
figure;
subplot(1,3,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
ax1 = scatter(T1.soma_center(:,1),T1.soma_center(:,2),3,color3,'filled');
set(gca, 'YDir','reverse');
xline(center(1),'--k');
yline(center(2),'--k');
axis image; axis off;
vector_scale = 2;
ring_scale = 10;
%%
range = [-pi,-pi/2; -pi/2,0; 0, pi/2; pi/2, pi];
%%
clear axon_terminals1 axon_points1 soma1
for kk = 1:4
    clear a1 axon_terminal_temp axon_points_temp soma_all_temp
    a1 = (T1.soma_angle>= range(kk,1) & T1.soma_angle<range(kk,2));
    axon_terminal_temp = axon_terminal_all3(a1);
    axon_points_temp = axon_points_all3(a1);
    soma_all_temp = soma_all3(a1);
    %%
    clear axon_terminals axon_points soma
    for i = 1:numel(soma_all_temp)
        clear axon_terminal_temp1 axon_points_temp1 soma_temp1
        axon_terminal_temp1 = double(axon_terminal_temp{i});
        axon_points_temp1 = double(axon_points_temp{i});
        soma_temp1 = double(soma_all_temp{i});
        
        axon_terminal_temp2 = zeros(size(axon_terminal_temp1));
        axon_points_temp2 = zeros(size(axon_points_temp1));
        axon_terminal_temp2 = axon_terminal_temp1-soma_temp1(1,:);
        axon_points_temp2 = axon_points_temp1-soma_temp1(1,:);
        axon_terminals{i} = axon_terminal_temp2;
        axon_points{i} = axon_points_temp2;
        soma{i} =  soma_temp1(1,:);
    end
    axon_terminals1{kk} = cat(1,axon_terminals{:});
    axon_points1{kk} = cat(1,axon_points{:});
    soma1{kk} = cat(1,soma{:});
end
%%
figure;
for kk = 1:4
    subplot(1,4,kk);
    clear axon_points1a axon_terminals1a
    axon_points1a = axon_points1{kk};
    axon_terminals1a = axon_terminals1{kk};
    soma1a = soma1{kk};
    soma_mean = mean(soma1a);
    overlayOutlines(coords,scale,[0.8,0.8,0.8]);
    set(gca, 'YDir','reverse');
    plot(axon_points1a(:,3)+soma_mean(3), axon_points1a(:,1)+soma_mean(1), '.', 'Color', 'k','MarkerSize',1); 
    hold on;
    plot(axon_terminals1a(:,3)+soma_mean(3), axon_terminals1a(:,1)+soma_mean(1), '.', 'Color', 'r','MarkerSize',5); 
%     xlim([-200,400]);
%     ylim([-300,400]);
    axis off; axis image;
    xlim([0,600]);
    ylim([0,1200]);
end
