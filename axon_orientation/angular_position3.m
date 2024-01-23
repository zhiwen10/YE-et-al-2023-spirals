%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
addpath(genpath(fullfile(githubDir, 'FMAToolbox')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\CircStat2012a'));
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
clear a1 T1
a1 = not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"])));
T1 = T(a1,:);
h = plot_distribution(T1,coords);
% print(h, 'SSp_axon_bias_version4', '-dpdf', '-bestfit', '-painters');
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
%%
pval = circ_wwtest(angle_diff_real,angle_diff_perm); % bad
%%
for i = 1:100
    [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end
%%
clear a1 T1
a1 = ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa"]));
T1 = T(a1,:);
h = plot_distribution(T1,coords);
print(h, 'peripheral_axon_bias_version4', '-dpdf', '-bestfit', '-painters');
%% calculate angular position of axon
scale = 1;
x1 = linspace(-pi,pi,180);
color2 = colorcet('C06','N',180);
color3 = interp1(x1,color2,T1.soma_angle);
center = [244,542]; % SSp-un
figure;
subplot(1,4,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
ax1 = scatter(T1.soma_center(:,1),T1.soma_center(:,2),3,color3,'filled');
set(gca, 'YDir','reverse');
xline(center(1),'--k');
yline(center(2),'--k');
axis image; axis off;
vector_scale = 2;
ring_scale = 10;

subplot(1,4,2);
r = 1;
z = r*exp(1i*T1.soma_angle);
ax1 = scatter(real(z)*ring_scale,imag(z)*ring_scale,3,color3,'filled');
set(gca, 'YDir','reverse');
axis image; axis off;
xlim([-20,20]);
ylim([-20,20]);
%
subplot(1,4,3);
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T1.bias_angle);

for icell = 1:numel(z)
    plot([real(z(icell))*ring_scale-T1.axon_bias(icell,1)*vector_scale*T1.pc_ratio(icell),...
        real(z(icell))*ring_scale+T1.axon_bias(icell,1)*vector_scale*T1.pc_ratio(icell)],...
        [imag(z(icell))*ring_scale-T1.axon_bias(icell,2)*vector_scale*T1.pc_ratio(icell),...
        imag(z(icell))*ring_scale+T1.axon_bias(icell,2)*vector_scale*T1.pc_ratio(icell)],...
        'color',colora3(icell,:),'LineWidth',1);
    hold on;
end
set(gca, 'YDir','reverse');
axis image; axis off;
xlim([-20,20]);
ylim([-20,20]);
%
subplot(1,4,4);
colormap(color2);
cb = colorbar;
axis off; 
colorbar('Ticks',[0,0.5,1],...
         'TickLabels',{'-\pi','0','\pi'})
%%
figure;
x1 = linspace(-90,90,180);
color2 = colorcet('C06','N',180);
color3 = interp1(x1,color2,T1.bias_angle);

edges1 = linspace(-pi,pi,5);
[N,edges1] = histcounts(T1.soma_angle,edges1);
for i = 1:numel(edges1)-1
    clear a ang_diff_temp z1 vector_all_sorted1 colora3_sorted1
    a = find(T1.soma_angle>=edges1(i) & T1.soma_angle<edges1(i+1));
    indx1{i} = a;
    % a = 1:numel(soma_angle_sorted);
    ang_diff_temp = T1.center_bias_angle(a);
    
    subplot(2,4,i);
    z1 = T1.soma_polar_angle(a);
    vector_all_sorted1 = T1.axon_bias(a,:);
    colora3_sorted1 = colora3(a,:);
    pc_ratio_sorted1 = T1.pc_ratio(a,:);
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
%%
figure; 
edges = [0:0.2:8];
histogram(T1.pc_ratio,edges);
%% circular bias is significant above pc_ratio>0.2
clear a1 T1
a1 = (not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"]))) & T.pc_ratio>=0.2 & T.pc_ratio<=0.5);
T1 = T(a1,:);
h = plot_distribution(T1,coords);
print(h, 'SSp_axon_bias_version4_ratio_lessthan_02_05', '-dpdf', '-bestfit', '-painters');
%%
clear a1 T1
a1 = (not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"]))) & T.pc_ratio>0.5);
T1 = T(a1,:);
h = plot_distribution(T1,coords);
print(h, 'SSp_axon_bias_version4_ratio_morethan_0.5', '-dpdf', '-bestfit', '-painters');
%% relationship of center distance vs pc_ratio
figure;
subplot(1,2,1)
clear orthog_vector soma_distance N pc_ratio_mean pc_ratio_sem
orthog_vector = T1.soma_center-center;
soma_distance = vecnorm(orthog_vector,2,2);
edges = [0:50:400];
[N,edges] = histcounts(soma_distance,edges);
for i = 1:numel(edges)-1
    clear a pc_ratio_temp
    a = find(soma_distance>=edges(i) & soma_distance<edges(i+1));
    %%
    pc_ratio_temp = T1.pc_ratio(a);
    pc_ratio_mean(i) = mean(pc_ratio_temp);
    pc_ratio_sem(i) = std(pc_ratio_temp)./sqrt(numel(pc_ratio_temp));
end
errorbar(edges(1:end-1),pc_ratio_mean,pc_ratio_sem);
%% mean circular_bias based on soma location
clear a1 T1
a1 = not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"])));
% a1 = (not(ismember(T.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"]))) & T.pc_ratio>=0.2);
T1 = T(a1,:);
%%
half_bin = pi/24;
bias_angle = deg2rad(T1.bias_angle);
bias_smooth_mean = zeros(numel(bias_angle),1);
bias_smooth_var = zeros(numel(bias_angle),1);
for kk = 1:size(T1,1)
    bias_temp = T1.soma_angle(kk);
    if bias_temp >= (pi-half_bin)
        min_range = bias_temp-half_bin;
        max_range = wrapToPi(bias_temp+half_bin);
        a = find(T1.soma_angle>=min_range | T1.soma_angle<=max_range);
    elseif bias_temp <= (-pi+half_bin)
        min_range = wrapToPi(bias_temp-half_bin);
        max_range = bias_temp+half_bin;
        a = find(T1.soma_angle>=min_range | T1.soma_angle<=max_range);
    else
        min_range = bias_temp-half_bin;
        max_range = bias_temp+half_bin;
        a = find(T1.soma_angle>=min_range & T1.soma_angle<max_range);
    end       
    soma_temp = T1.soma_angle(a);
    weight = T1.pc_ratio(a);
    bias_angle_temp = bias_angle(a);
    bias_angle_mean1 = circ_mean(bias_angle_temp);
    bias_angle_var1 =  circ_var(bias_angle_temp);
    bias_smooth_mean(kk) = bias_angle_mean1;
    bias_smooth_var(kk) = bias_angle_var1;
end
    %%
% for i = 1:numel(edges)-1
%     clear a bias_angle_mean1 bias_angle_temp
%     a = find(T1.soma_angle>=edges(i) & T1.soma_angle<edges(i+1));
%     bias_angle_temp = bias_angle(a);
%     bias_angle_mean1 = circ_mean(bias_angle_temp);
%     bias_angle_mean(i) = bias_angle_mean1;
%     bias_smooth(a) = bias_angle_mean1;
% end
%%
x1 = linspace(-pi,pi,180);
color2 = colorcet('C06','N',180);
color3 = interp1(x1,color2,T1.soma_angle);

x2 = linspace(-pi/2,pi/2,180);
color4 = colorcet('C06','N',180);
color5= interp1(x2,color4,bias_smooth_mean);

x3 = linspace(0,1,180);
color6 = parula(180);
color7= interp1(x3,color6,bias_smooth_var);

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
subplot(1,3,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
ax1 = scatter(T1.soma_center(:,1),T1.soma_center(:,2),3,color5,'filled');
set(gca, 'YDir','reverse');
xline(center(1),'--k');
yline(center(2),'--k');
axis image; axis off;
subplot(1,3,3);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
ax1 = scatter(T1.soma_center(:,1),T1.soma_center(:,2),3,color7,'filled');
set(gca, 'YDir','reverse');
xline(center(1),'--k');
yline(center(2),'--k');
axis image; axis off;