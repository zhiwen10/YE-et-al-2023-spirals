%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
ctr=0; 
% for nIdx = 1:numel(allCoords)
for nIdx =50
    somaCoords = allCoords{nIdx}(:,1)==1;
    coordsIdx = round(allCoords{nIdx}(:,2:4)/10); % compatibility with av (10 micron resolution)
    soma = [coordsIdx(somaCoords,1),  coordsIdx(somaCoords,2),  coordsIdx(somaCoords,3)];
    avIdxSoma = sub2ind(size(av), soma(1),soma(2),soma(3));
    stIdxSoma = av(avIdxSoma);
end
%%
figure;
for i = 1
    ax1 = subplot(1,1,1)
    % ax1 = subplot(2,5,i)
    plotBrainGrid([],ax1);
    index = i;
    noi = allCoords{index}; % neuron of interest: brain ID 210254, neuron ID *20600
    hold on;
    colors = get(gca, 'ColorOrder');
    for type = 1:max(noi(:,1))
        axonCoords = noi(:,1)==type;

        plot3(noi(axonCoords,2)/10, noi(axonCoords,4)/10, noi(axonCoords,3)/10, '.', 'Color', colors(type,:)); % y and z coords need to be flipped around here, I think this may have to do with plot3 and BrainGrid compatibility?
        hold on;
    end
    axis equal
    text(0,0,num2str(index))
    view(2)
end
%%
region_label = {"SSp-ll"};
st_region_indx = get_region_index_list(region_label,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
%%
axon_region_label = {"Isocortex"};
axon_st_cortex_indx = get_region_index_list(axon_region_label,st); % get all index values that contain "Isocortex" in st table
%%
for i = 10
    figure;
    ax1 = subplot(1,1,1);
    plotBrainGrid([],ax1);
    index = i;
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    hold on;
    colors = get(gca, 'ColorOrder');
    % for type = 1:max(noi(:,1))
    type = 2;
    axonCoords = noi(noi(:,1)==type,:);
    % filter axon terminal in cortex
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx);
    plot3(axonCoords_filt(:,1), axonCoords_filt(:,3), axonCoords_filt(:,2), '.', 'Color', colors(type,:)); 
    axis equal
    text(0,0,num2str(index))
    view(2)
end

%%
type = 2;
axon_all1 = [];
for i = 1:10
    clear noi axonCoords axonCoords_filt
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    axonCoords = noi(noi(:,1)==type,:);
    % filter axon terminal in cortex
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx);
    axon_all1 = [axon_all1;axonCoords_filt];
end

%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);

%% plot all axons in a region
regions_all = {"SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","VISp"};
cell_limit = 5;
for i = 1:7
    % region_label = {"SSp-ll"};
    region_label = regions_all{i};
    % [soma_all{i},axon_all{i}] = get_all_axon_in_area(st,av,allCoords,region_label,axon_st_cortex_indx,cell_limit);
    [soma_all{i},axon_all{i}] = get_all_axon_in_area(st,av,allCoords,region_label,axon_st_cortex_indx);
    % [soma_all{i},axon_all{i}] = get_all_axon_in_area2(st,av,allCoords,region_label);
end
color1 = hsv(7);

scale = 1;
figure;
ax1 = subplot(2,4,1)
overlayOutlines(coords,scale);
set(gca, 'YDir','reverse');
hold on;
for i = 1:6
    soma_current = soma_all{i};
    plot(ax1,soma_current(:,3), soma_current(:,1), '.', 'Color', color1(i,:),'MarkerSize',12); 
    axis off; axis image;
end

for i = 1:6 
    clear axon_current
    axon_current = axon_all{i};
    axon_current_2d = double(axon_current(:,[3,1]));
    center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(center,size(axon_current_2d,1),1);
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    
    subplot(2,4,i+1)
    overlayOutlines(coords,scale);
    set(gca, 'YDir','reverse');
    hold on;
    plot(axon_current(:,3), axon_current(:,1), '.', 'Color', color1(i,:),'MarkerSize',2); 
    hold on;
    plot([center(1)-U(1,1)*50,center(1)+U(1,1)*50],[center(2)-U(2,1)*50,center(2)+U(2,1)*50],'k','LineWidth',4)
    axis off; axis image;
    title(regions_all{i});
    
    ax2 = subplot(2,4,8);
    overlayOutlines(coords,scale);
    set(gca, 'YDir','reverse');
    hold on;
    scatter(ax2,center(1),center(2),20,'MarkerFaceColor',color1(i,:))
    plot(ax2,[center(1)-U(1,1)*50,center(1)+U(1,1)*50],[center(2)-U(2,1)*50,center(2)+U(2,1)*50],'color',color1(i,:),'LineWidth',4)
    axis off; axis image;
    
    
end
%%
coords_current = axon_all{4};
axon_current_2d = double(coords_current(:,[3,1]));
center = mean(axon_current_2d,1);
axon_current_2d = axon_current_2d-repmat(center,size(axon_current_2d,1),1);
[U,S,V] = svd(double(axon_current_2d)',"econ");
figure;
scatter(center(1)+axon_current_2d(:,1),center(2)+axon_current_2d(:,2))
hold on;
plot([center(1)-U(1,1)*100,center(1)+U(1,1)*100],[center(2)-U(1,2)*100,center(2)+U(1,2)*100],'r')
hold on;
plot([center(1)-U(2,1)*100,center(1)+U(2,1)*100],[center(2)-U(2,2)*100,center(2)+U(2,2)*100],'k')
axis equal
set(gca, 'YDir','reverse')
%%
figure;
indx = randperm(size(axon_current_2d,1),10000);
rand_2d = axon_current_2d(indx,:);
scatter_kde(center(1)+rand_2d(:,1),center(2)+rand_2d(:,2), 'filled');
hold on;
plot([center(1),center(1)+U(1,1)*100],[center(2),center(2)+U(1,2)*100],'r')
hold on;
plot([center(1),center(1)+U(2,1)*100],[center(2),center(2)+U(2,2)*100],'k')
axis equal
set(gca, 'YDir','reverse')
%%
region_label = {"CLA"};
st_region_indx = get_region_index_list(region_label,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
%%
color1 = hsv(17);
figure;
ax1 = subplot(1,2,1);
plotBrainGrid([],ax1);
for i = 1:17
    clear noi axonCoords axonCoords_filt
    index = i;
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    hold on;
    type = 1;
    somaCoords = noi(noi(:,1)==type,2:4);
    plot3(somaCoords(:,1)/10, somaCoords(:,3)/10, somaCoords(:,2)/10, '.', 'Color', color1(i,:),'MarkerSize',12); 
    axis equal
end
view(2)
title('soma')

ax1 = subplot(1,2,2);
plotBrainGrid([],ax1);
for i = 1:17
    clear noi axonCoords axonCoords_filt
    index = i;
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    hold on;
    type = 2;
    axonCoords = noi(noi(:,1)==type,:);
    % filter axon terminal in cortex
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx);
    plot3(axonCoords_filt(:,1), axonCoords_filt(:,3), axonCoords_filt(:,2), '.', 'Color', color1(i,:),'MarkerSize',2); 
    axis equal
end
view(2)
title('axon')
%%
f = figure; 
for i = 1:17
    ax1 = subplottight(4,5,i)
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    plotBrainGrid([],ax1); hold on;
    colors = get(gca, 'ColorOrder');
    for type = 1:max(noi(:,1))
        theseCoords = noi(:,1)==type;
        if type ==1
            markerS = 12;
        else
            markerS = 2;
        end
        plot3(noi(theseCoords,2)/10, noi(theseCoords,4)/10, noi(theseCoords,3)/10, '.', 'Color', colors(type,:),'MarkerSize',markerS); % y and z coords need to be flipped around here, I think this may have to do with plot3 and BrainGrid compatibility?
        hold on;
    end
    axis equal
%     view(3)
end