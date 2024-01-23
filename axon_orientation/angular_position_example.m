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
%% cortex surface outline
root1 = {'/997/'};
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
cell_id = 40;
% for k = 1:numel(regions_all1)
for k = 1
    current_axon = axon_points_all2{k};
    current_soma = soma_all2{k};
    current_terminal = axon_terminal_all2{k};
    cell_n = numel(current_axon);
    % only select cells that have axon branching points
    a = cellfun(@(x) not(isempty(x)), current_terminal);
    current_terminal = current_terminal(a);
    current_soma = current_soma(a);

    [soma_center,vector_all, angle1, polarity] = get_axon_pc(current_terminal,current_soma);   
    color2= colorcet('C06','N',180);
    angle1(angle1==91) = 90;
    color3 = color2(angle1+90,:);
%     if cell_n>40
%         cell_n =40;
%     end
    
    h = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
    % for icell = 1:cell_n
    for icell = 1:7
        clear axon_current_all
        icella = cell_id+icell;
        current_axon_all = current_axon{icella}/5;
        current_terminal_all = current_terminal{icella}/5;
        current_soma_all = current_soma{icella}/5;       
        subplot(8,5,icell);
        scale3 = 1;
        plotOutline_fill(root1,st,atlas1,[],scale3,'w');
        axis image; axis off;
        set(gca, 'YDir','reverse');              
        hold on;
        scatter(current_axon_all(:,3), current_axon_all(:,1),0.5,'MarkerFaceColor','k','MarkerEdgeColor','None'); 
        hold on;
        scatter(current_soma_all(:,3),current_soma_all(:,1),2,'r','filled');
            hold on;
        plot([soma_center(icella,1)/5-vector_all(icella,1)*20*polarity(icella),soma_center(icella,1)/5+vector_all(icella,1)*20*polarity(icella)],...
            [soma_center(icella,2)/5-vector_all(icella,2)*20*polarity(icella),soma_center(icella,2)/5+vector_all(icella,2)*20*polarity(icella)],...
            'color',color3(icella,:),'LineWidth',1)
    
        axis off; axis image;
    end
    print(h, [char(regions_all1{k}) '-' num2str(cell_id(1)+1)], '-dpdf', '-bestfit', '-painters');
end
%% SSp
area = [8,9,10,11,12,13]; % 'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd'
example_id = [1,6,2,32,1,52];
h1 = figure('Renderer', 'painters', 'Position', [100 100 700 500]);
subplot(2,2,1);
scale3 = 1; scale=5;
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
    scatter(current_axon{icell}(:,3)/5, current_axon{icell}(:,1)/5,0.5,'MarkerFaceColor',colora3(icell,:),'MarkerEdgeColor','None'); 
    hold on;
    plot([soma_center(icell,1)/5-vector_all(icell,1)*20*polarity(icell),soma_center(icell,1)/5+vector_all(icell,1)*20*polarity(icell)],...
        [soma_center(icell,2)/5-vector_all(icell,2)*20*polarity(icell),soma_center(icell,2)/5+vector_all(icell,2)*20*polarity(icell)],...
        'color','k','LineWidth',1)
    axis off; axis image;
end
%%
clear area example_id
% area = [8,9,10,11,12,13,15,3,1,2]; %'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd', 'SSs','AUD','TEa','VIS','RSP'
% example_id = [1,6,2,32,1,52,39,18,2,15];
area = [15,3,1,2]; %'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd', 'SSs','AUD','TEa','VIS','RSP'
example_id = [39,18,2,15];
subplot(2,2,2);
hold off;
scale3 = 1; scale=5;
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
    scatter(current_axon{icell}(:,3)/5, current_axon{icell}(:,1)/5,0.5,'MarkerFaceColor',colora3(icell,:),'MarkerEdgeColor','None'); 
    hold on;
    plot([soma_center(icell,1)/5-vector_all(icell,1)*20*polarity(icell),soma_center(icell,1)/5+vector_all(icell,1)*20*polarity(icell)],...
        [soma_center(icell,2)/5-vector_all(icell,2)*20*polarity(icell),soma_center(icell,2)/5+vector_all(icell,2)*20*polarity(icell)],...
        'color','k','LineWidth',1)
    axis off; axis image;
end
