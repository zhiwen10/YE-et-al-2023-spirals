%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
load('all_cell_with_parents.mat');
%%
axon_region_label = {"MOs","MOp","ACA"};
axon_st_cortex_indx = get_region_index_list(axon_region_label,st); % get all index values that contain "Isocortex" in st table

axon_region_label2 = {"CTX","cc"};
axon_st_cortex_indx2 = get_region_index_list(axon_region_label2,st); % get all index values that contain "Isocortex" in st table
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
regions_all1 = {"MOs","MOp","ACA"};
hemi_terminal = 1;
hemi_axon_all = 2;
axon_points_all2 = {};
soma_all2 = {};
axon_terminal_all2 = {};
for k = 1:numel(regions_all1)
% for k = 2
    iregion = regions_all1{k};
    st_region_indx = get_region_index_list(iregion,st); % get all index values that contain "SSp-ll" in st table
    [cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
    % only get terminals in the left sensory cortex
    [~,axon_terminal_all1] = get_axon_branching_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx, hemi_terminal);
    % get axon points for both hemispheres in the entire cortex 
    [soma_all1,axon_points_all1,dendrite_all1] = get_axon_points_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx, hemi_terminal);
    axon_points_all2{k} = axon_points_all1;
    soma_all2{k} = soma_all1;
    dendrite_all2{k} = dendrite_all1;
    axon_terminal_all2{k} = axon_terminal_all1;
end
%%
cell_n = cellfun(@(x) numel(x),axon_points_all2);
%%
color4 = hsv(numel(regions_all1));
for k = 1:numel(regions_all1)
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
    if cell_n>40
        cell_n =40;
    end
    
    h = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
    for icell = 1:cell_n
        clear axon_current_all
        current_axon_all = current_axon{icell}/5;
        current_terminal_all = current_terminal{icell}/5;
        current_soma_all = current_soma{icell}/5;       
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
        plot([soma_center(icell,1)/5-vector_all(icell,1)*20*polarity(icell),soma_center(icell,1)/5+vector_all(icell,1)*20*polarity(icell)],...
            [soma_center(icell,2)/5-vector_all(icell,2)*20*polarity(icell),soma_center(icell,2)/5+vector_all(icell,2)*20*polarity(icell)],...
            'color',color3(icell,:),'LineWidth',1)
    
        axis off; axis image;
    end
    savefig(h,[char(regions_all1{k}) '-' num2str(cell_id(1))]);
    print(h, [char(regions_all1{k}) '-' num2str(cell_id(1))], '-dpdf', '-bestfit', '-painters');
end
%%
iregion = regions_all1{2};
st_region_indx = get_region_index_list(iregion,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
icell = 8;
noi = allCoords{cell_id(icell)};
type = 2;
axonCoords1 = noi(noi(:,1)== type,:);  
plotBrainGrid; hold on;
plot3(axonCoords1(:,2)/10, axonCoords1(:,4)/10, axonCoords1(:,3)/10, '.', 'Color', 'r'); 
% plot3(current_axon_all(:,1), current_axon_all(:,3), current_axon_all(:,2), '.', 'Color', 'r'); 
axis equal
%%
current_axon = axon_points_all2{2};
current_axon_all = current_axon{1};
plotBrainGrid; hold on;
plot3(current_axon_all(:,1), current_axon_all(:,3), current_axon_all(:,2), '.', 'Color', 'r'); 
axis equal
%% plot all terminals by area
regions_all1 = {"MOs","MOp","ACA"};
maskPath{1} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{2} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
%
id = [1,2,3];
color4 = hsv(14);
hemi = [];
hemi1 = 1;
h = figure('Renderer', 'painters', 'Position', [100 100 500 500]);
count1 = 1;
for k=1:3
    idi = id(k);
    current_axon = axon_points_all2{idi};
    current_soma = soma_all2{idi};
    current_dendrite = dendrite_all2{idi};
    current_terminal = axon_terminal_all2{idi};
    soma_n = numel(current_soma);
    soma1 = cellfun(@(x) x(1,:),current_soma,'UniformOutput',false);
    soma2 = cat(1,soma1{:});
    [a,b] = sort(soma2(:,1));
    ax{k} = subplottight(3,3,count1);
    scale3 = 1;
    plotOutline_fill(root1,st,atlas1,hemi,scale3,'w');
    for kk = 1:3
        plotOutline_fill(maskPath(id(kk)),st,atlas1,hemi,scale3,'w');
    end
    axis image; axis off;
    current_axon_all = cat(1,current_axon{:})/5;
    current_soma_all = cat(1,current_soma{:})/5;      
    current_terminal_all = cat(1,current_terminal{:})/5;  
    scatter_kde(current_terminal_all(:,3), current_terminal_all(:,1),'filled','MarkerSize',0.5);
    
    faceAlpha =0.2;
    plotOutline_fill(maskPath(idi),st,atlas1,hemi1,scale3,[0.2,0.2,0.2],faceAlpha);
    count1 = count1+1;
    set(gca, 'YDir','reverse');
    xlim([-20,600]);  
    text(100,-50,regions_all1{idi});
end
%%
st_region_indx = get_region_index_list(regions_all1,st); 
[cell_id] = get_cell_id(allCoords,av,st_region_indx); 
% only get terminals in the left sensory cortex
[~,axon_terminal_all1] = get_axon_branching_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx, hemi_terminal);
% get axon points for both hemispheres in the entire cortex 
[soma_all1,axon_points_all1,dendrite_all1] = get_axon_points_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx, hemi_terminal);

a = cellfun(@(x) not(isempty(x)), axon_terminal_all1);
all_terminal = axon_terminal_all1(a);
all_soma = soma_all1(a);
[soma_center,vector_all, angle1, polarity] = get_axon_pc(all_terminal,all_soma);   
%%
center = [244,542]; % SSp-un
orthog_vector = soma_center-center;
for i =1:size(vector_all,1)
    u = [vector_all(i,:) 0];
    v = [orthog_vector(i,:) 0];
    ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
end    
ang_diff2 = rad2deg(ang_diff1);
edges = 0:10:180;
[N1,edges] = histcounts(ang_diff2,edges);
h1 = figure;
subplot(2,1,1);
stairs(edges(1:end-1),N1,'r');
subplot(2,1,2)
scatter(ang_diff2,polarity,'k','filled')
%%
id = [1,2,3];
count1 = 1;
h1 = figure;
for k = 1:numel(id)
    clear soma_center vector_all angle1 polarity orthog_vector u v ang_diff1 ang_diff2 N1
    current_soma = soma_all2{id(k)};
    current_terminal = axon_terminal_all2{id(k)};
    
    a = cellfun(@(x) not(isempty(x)), current_terminal);
    current_terminal = current_terminal(a);
    current_soma = current_soma(a);

    [soma_center,vector_all, angle1, polarity] = get_axon_pc(current_terminal,current_soma);   

    orthog_vector = soma_center-center;
    for i =1:size(vector_all,1)
        u = [vector_all(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N1,edges] = histcounts(ang_diff2,edges);

    subplot(3,2,2*k-1);
    stairs(edges(1:end-1),N1,'r');
    title(regions_all1{id(k)});
    subplot(3,2,2*k)
    scatter(ang_diff2,polarity,'k','filled')
    count1 = count1+1;
end