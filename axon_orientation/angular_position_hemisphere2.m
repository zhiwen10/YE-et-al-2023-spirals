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
axon_region_anterior = {'MOp','MOs','ACA'};
anterior_indx= get_region_index_list(axon_region_anterior,st); % get all index values that contain "Isocortex" in st table
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
hemi_terminal = 2;
hemi_axon_all = [];
axon_points_all2 = {};
soma_all2 = {};
axon_terminal_all2 = {};
for k = 1:numel(regions_all1)
% for k = 2
    iregion = regions_all1(k);
    st_region_indx = get_region_index_list(iregion,st); % get all index values that contain "SSp-ll" in st table
    [cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
    % only get terminals in the left sensory cortex
    [~,axon_terminal_all1] = get_axon_branching_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx2, hemi_terminal);
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
axon_points_all3 = axon_points_all2;
% soma_all3 = soma_all2([8:14]);
% axon_terminal_all3 = axon_terminal_all2([8:14]);
soma_all3 = cat(2, soma_all3{:});
axon_terminal_all3 = cat(2,axon_terminal_all3{:});
axon_points_all3 = cat(2,axon_points_all3{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
axon_points_all3 = axon_points_all3(a);
%%
cell_id = 1;
for k = 1:numel(regions_all1)
    current_axon = axon_points_all2{k};
    current_soma = soma_all2{k};
    current_terminal = axon_terminal_all2{k};
    % only select cells that have axon branching points
    a = cellfun(@(x) not(isempty(x)), current_terminal);
    current_terminal = current_terminal(a);
    current_soma = current_soma(a);
    cell_n = numel(current_soma);

    [soma_center,vector_all, angle1, polarity] = get_axon_pc(current_terminal,current_soma);   
    color2= colorcet('C06','N',180);
    angle1(angle1==91) = 90;
    color3 = color2(angle1+90,:);
    if cell_n>40
        cell_n =40;
    end
    
    h = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
    for icell = 1:cell_n
        clear axon_current_all
        current_axon_all = current_axon{icell}/5;
        current_terminal_all = current_terminal{icell}/5;
        current_soma_all = current_soma{icell}/5;   
        axon_center = mean(current_terminal_all,1);
        subplot(8,5,icell);
        scale3 = 1;
        plotOutline_fill(root1,st,atlas1,[],scale3,'w');
        axis image; axis off;
        set(gca, 'YDir','reverse');
%         xlim([-20,axis_lim(1,iplot)]+20);
%         ylim([-20,axis_lim(2,iplot)]+20);                
        hold on;
        scatter(current_axon_all(:,3), current_axon_all(:,1),0.5,'MarkerFaceColor','k','MarkerEdgeColor','None'); 
%         hold on;
%         scatter(current_terminal_all(:,3),current_terminal_all(:,1),0.5,'MarkerFaceColor','k','MarkerEdgeColor','None'); 
        hold on;
        scatter(current_soma_all(:,3),current_soma_all(:,1),2,'r','filled');
            hold on;
        plot([axon_center(3)-vector_all(icell,1)*5*polarity(icell),axon_center(3)+vector_all(icell,1)*5*polarity(icell)],...
            [axon_center(1)-vector_all(icell,2)*5*polarity(icell),axon_center(1)+vector_all(icell,2)*5*polarity(icell)],...
            'color',color3(icell,:),'LineWidth',1)
    
        axis off; axis image;
    end
    savefig(h,[char(regions_all1{k}) '-' num2str(cell_id(1))]);
    print(h, [char(regions_all1{k}) '-' num2str(cell_id(1))], '-dpdf', '-bestfit', '-painters');
end

%%
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
T = get_axon_bias_table(soma_center,vector_all,angle1,polarity);
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = T.soma_center;
cell_n = size(T,1);
%%
scale = 1;
color2= colorcet('C06','N',180);
x1 = linspace(-pi,pi,180);
color3 = interp1(x1,color2,T.soma_angle);
h1 = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,1200]);
ylim([0,1200]);
for k = 1:numel(axon_terminal_all3)
    current_axon = axon_points_all3{k};
    current_terminal = axon_terminal_all3{k};
    hold on;
    scatter(soma_center1(k,1),soma_center1(k,2),6,'MarkerFaceColor',color3(k,:),'MarkerEdgeColor','None'); 
    hold on;
    scatter(current_terminal(:,3)/scale, current_terminal(:,1)/scale,1,'MarkerFaceColor',color3(k,:),'MarkerEdgeColor','None'); 
    hold on;
end
%%
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
T = get_axon_bias_table(soma_center,vector_all,angle1,polarity);
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = T.soma_center;
cell_n = size(T,1);

h1 = figure('Renderer', 'painters', 'Position', [100 100 700 500]);
ax1 = subplot(2,3,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([500,1100]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-T.axon_bias(i,1)*scale1*T.pc_ratio(i),soma_center1(i,1)+T.axon_bias(i,1)*scale1*T.pc_ratio(i)],...
        [soma_center1(i,2)-T.axon_bias(i,2)*scale1*T.pc_ratio(i),soma_center1(i,2)+T.axon_bias(i,2)*scale1*T.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.05 a(2)+0.1 a(3) a(4)/2])% To change size

ax2 = subplot(2,3,3);
edges = 0:10:180;
[N1,edges] = histcounts(T.center_bias_angle,edges);
% stairs(edges(1:end-1),N1,'r');
% xlim([0,180]);
% ylim([0,40]);

% shuffle
clear ang_diff1 ang_diff2
center = [244,542]; % SSp-un
orthog_vector = soma_center1-center;
for k = 1:100
    rp_indx = randperm(size(T.axon_bias,1));
    vector_all1 = T.axon_bias(rp_indx,:);
    for i =1:size(T.axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
ax2 = plot_angle_stairs_shuffle(ax2,edges,N1,N_rp);
xlim([0,180]); ylim([0,40]);
%%
soma_all3 = soma_all2;
axon_terminal_all3 = axon_terminal_all2;
soma_all3 = cat(2, soma_all3{:});
axon_terminal_all3 = cat(2,axon_terminal_all3{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
[soma_center,vector_all, angle1, polarity] = get_axon_pc(axon_terminal_all3,soma_all3);   
T = get_axon_bias_table(soma_center,vector_all,angle1,polarity);
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T.bias_angle);

vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = T.soma_center;
cell_n = size(T,1);
ax3 = subplot(2,3,5);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([500,1100]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-T.axon_bias(i,1)*scale1*T.pc_ratio(i),soma_center1(i,1)+T.axon_bias(i,1)*scale1*T.pc_ratio(i)],...
        [soma_center1(i,2)-T.axon_bias(i,2)*scale1*T.pc_ratio(i),soma_center1(i,2)+T.axon_bias(i,2)*scale1*T.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.05 a(2)+0.1 a(3) a(4)/2])% To change size

ax4 = subplot(2,3,6);
edges = 0:10:180;
[N1,edges] = histcounts(T.center_bias_angle,edges);   
% shuffle
center = [244,542]; % SSp-un
orthog_vector = soma_center1-center;
for k = 1:100
    rp_indx = randperm(size(T.axon_bias,1));
    vector_all1 = T.axon_bias(rp_indx,:);
    for i =1:size(T.axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
ax4 = plot_angle_stairs_shuffle(ax4,edges,N1,N_rp);
xlim([0,180]); ylim([0,40]);
%% SSp
area = [8,9,10,11,11,12,13,13,13]; % 'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd'
example_id = [1,6,2,32,39,20,52,39,21];
subplot(2,3,1);
hold off;
scale3 = 1; scale=1;
scale1 = 15;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
axis image; axis off;
set(gca, 'YDir','reverse'); 
for areai = 1:numel(area)
    k1 = area(areai);
    icell = example_id(areai);
    current_axon = axon_points_all2{k1};
    current_soma = soma_all2{k1};
    current_terminal = axon_terminal_all2{k1};
    cell_n = numel(current_axon);    
    a = cellfun(@(x) not(isempty(x)), current_terminal);
    current_terminal = current_terminal(a);
    current_soma = current_soma(a);
    
    [soma_center,vector_all, angle1, polarity] = get_axon_pc(current_terminal,current_soma);   
    T = get_axon_bias_table(soma_center,vector_all,angle1,polarity);
    colora2= colorcet('C06','N',180);
    x1 = linspace(-90,90,180);
    colora3 = interp1(x1,colora2,T.bias_angle);    
    hold on;
    scatter(current_axon{icell}(:,3)/scale, current_axon{icell}(:,1)/scale,0.5,'MarkerFaceColor',colora3(icell,:),'MarkerEdgeColor','None'); 
    hold on;
    plot([soma_center(icell,1)/scale-vector_all(icell,1)*scale1*polarity(icell),soma_center(icell,1)/scale+vector_all(icell,1)*scale1*polarity(icell)],...
        [soma_center(icell,2)/scale-vector_all(icell,2)*scale1*polarity(icell),soma_center(icell,2)/scale+vector_all(icell,2)*scale1*polarity(icell)],...
        'color','k','LineWidth',1)
    axis off; axis image;
end
xlim([0,600]);
ylim([0,1200]);
%
clear area example_id
% area = [8,9,10,11,12,13,15,3,1,2]; %'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd', 'SSs','AUD','TEa','VIS','RSP'
% example_id = [1,6,2,32,1,52,39,18,2,15];
area = [15,3,1,2]; % 'SSs','AUD','TEa','VIS','RSP'
example_id = [27,18,2,15];
subplot(2,3,4);
hold off;
scale3 = 1;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
axis image; axis off;
set(gca, 'YDir','reverse'); 
for areai = 1:numel(area)
    k1 = area(areai);
    icell = example_id(areai);
    current_axon = axon_points_all2{k1};
    current_soma = soma_all2{k1};
    current_terminal = axon_terminal_all2{k1};
    cell_n = numel(current_axon);    
    a = cellfun(@(x) not(isempty(x)), current_terminal);
    current_terminal = current_terminal(a);
    current_soma = current_soma(a);
    
    [soma_center,vector_all, angle1, polarity] = get_axon_pc(current_terminal,current_soma);   
    T = get_axon_bias_table(soma_center,vector_all,angle1,polarity);
    colora2= colorcet('C06','N',180);
    x1 = linspace(-90,90,180);
    colora3 = interp1(x1,colora2,T.bias_angle);    
    hold on;
    scatter(current_axon{icell}(:,3)/scale, current_axon{icell}(:,1)/scale,0.5,'MarkerFaceColor',colora3(icell,:),'MarkerEdgeColor','None'); 
    hold on;
    plot([soma_center(icell,1)/scale-vector_all(icell,1)*scale1*polarity(icell),soma_center(icell,1)/scale+vector_all(icell,1)*scale1*polarity(icell)],...
        [soma_center(icell,2)/scale-vector_all(icell,2)*scale1*polarity(icell),soma_center(icell,2)/scale+vector_all(icell,2)*scale1*polarity(icell)],...
        'color','k','LineWidth',1)
    axis off; axis image;
end
xlim([0,600]);
ylim([0,1200]);
%%
print(h1, 'axon_bias_all', '-dpdf', '-bestfit', '-painters');