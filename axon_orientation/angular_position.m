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
end
%%
soma_all3 = soma_all2;
axon_terminal_all3 = axon_terminal_all2;

soma_all3 = cat(2, soma_all3{:});
axon_terminal_all3 = cat(2,axon_terminal_all3{:});

a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
%%
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
%% calculate angular position of soma
clear ang_diff1 ang_diff2
center = [244,542]; % SSp-un
orthog_vector = soma_center-center;
for i =1:size(vector_all,1)
    u = [vector_all(i,:) 0];
    v = [orthog_vector(i,:) 0];
    ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
end    
ang_diff2 = rad2deg(ang_diff1);
%% calculate angular position of axon
vector_scale = 2;
ring_scale = 10;
soma_angle = atan2(orthog_vector(:,2),orthog_vector(:,1));
x1 = linspace(-pi,pi,180);
color2 = colorcet('C06','N',180);
color3 = interp1(x1,color2,soma_angle);
figure;
subplot(1,4,1);
ax1 = scatter(orthog_vector(:,1),orthog_vector(:,2),12,color3,'filled');
set(gca, 'YDir','reverse');
xline(0,'--k');
yline(0,'--k');
axis image;

subplot(1,4,2);
r = 1;
z = r*exp(1i*soma_angle);
ax1 = scatter(real(z)*ring_scale,imag(z)*ring_scale,12,color3,'filled');
set(gca, 'YDir','reverse');
axis image; axis off;
xlim([-20,20]);
ylim([-20,20]);
%
subplot(1,4,3);
colora2= colorcet('C06','N',180);
angle1(angle1==91) = 90;
colora3 = colora2(angle1+90,:);
for icell = 1:numel(z)
    plot([real(z(icell))*ring_scale-vector_all(icell,1)*vector_scale*polarity(icell),...
        real(z(icell))*ring_scale+vector_all(icell,1)*vector_scale*polarity(icell)],...
        [imag(z(icell))*ring_scale-vector_all(icell,2)*vector_scale*polarity(icell),...
        imag(z(icell))*ring_scale+vector_all(icell,2)*vector_scale*polarity(icell)],...
        'color',colora3(icell,:),'LineWidth',1);
    hold on;
end
set(gca, 'YDir','reverse');
axis image; axis off;
xlim([-20,20]);
ylim([-20,20]);

subplot(1,4,4);
colormap(color2);
cb = colorbar;
axis off; 
colorbar('Ticks',[0,0.5,1],...
         'TickLabels',{'-\pi','0','\pi'})
%% sort neuron by soma position
T = table(soma_center,soma_angle,z,vector_all,angle1,polarity,ang_diff2);
T = renamevars(T,["z","vector_all","angle1","polarity","ang_diff2"], ...
                 ["soma_polar_angle","axon_bias","bias_angle","pc_ratio","center_bias_angle"]);
%%
writetable(T,"axon_bias_table.csv");
%%
T1 = sortrows(T,'soma_angle');
% T1 = T;
soma_angle_sorted = T1.soma_angle;
soma_polar_angle_sorted = T1.soma_polar_angle;
bias_angle_sorted = T1.bias_angle;
axon_bias_sorted = T1.axon_bias;
center_bias_angle_sorted = T1.center_bias_angle;
pc_ratio_sorted = T1.pc_ratio;
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,bias_angle_sorted);
%%
soma_center1 = T1.soma_center;
cell_n = size(T1,1);
scale = 1;
figure;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias_sorted(i,1)*scale1*pc_ratio_sorted(i),soma_center1(i,1)+axon_bias_sorted(i,1)*scale1*pc_ratio_sorted(i)],...
        [soma_center1(i,2)-axon_bias_sorted(i,2)*scale1*pc_ratio_sorted(i),soma_center1(i,2)+axon_bias_sorted(i,2)*scale1*pc_ratio_sorted(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
% scatter(soma_center(:,1),soma_center(:,2),[],color3,'filled');
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
title('All Layers');
%%
figure;
edges = 0:10:180;
[N1,edges] = histcounts(ang_diff2,edges);
stairs(edges(1:end-1),N1,'r');
%%
figure;
edges1 = linspace(-pi,pi,5);
[N,edges1] = histcounts(soma_angle_sorted,edges1);
for i = 1:numel(edges1)-1
    clear a ang_diff_temp z1 vector_all_sorted1 colora3_sorted1
    a = find(soma_angle_sorted>=edges1(i) & soma_angle_sorted<edges1(i+1));
    indx1{i} = a;
    % a = 1:numel(soma_angle_sorted);
    ang_diff_temp = center_bias_angle_sorted(a);
    
    subplot(2,4,i);
    z1 = soma_polar_angle_sorted(a);
    vector_all_sorted1 = axon_bias_sorted(a,:);
    colora3_sorted1 = colora3(a,:);
    pc_ratio_sorted1 = pc_ratio_sorted(a,:);
    ax1 = scatter(real(z)*ring_scale,imag(z)*ring_scale,3,'k','filled');
    set(gca, 'YDir','reverse');
    hold on;
    for icell = 1:numel(z1)
        plot([real(z1(icell))*ring_scale-vector_all_sorted1(icell,1)*vector_scale*pc_ratio_sorted1(icell),...
            real(z1(icell))*ring_scale+vector_all_sorted1(icell,1)*vector_scale*pc_ratio_sorted1(icell)],...
            [imag(z1(icell))*ring_scale-vector_all_sorted1(icell,2)*vector_scale*pc_ratio_sorted1(icell),...
            imag(z1(icell))*ring_scale+vector_all_sorted1(icell,2)*vector_scale*pc_ratio_sorted1(icell)],...
            'color',colora3_sorted1(icell,:),'LineWidth',1);
        hold on;
    end
    set(gca, 'YDir','reverse');
    axis image; axis off;
    xlim([-20,20]);
    ylim([-20,20]);
    
    subplot(2,4,i+4);
    edges = 0:10:180;
    [N1,edges] = histcounts(ang_diff_temp,edges);
    stairs(edges(1:end-1),N1,'r');
end
    