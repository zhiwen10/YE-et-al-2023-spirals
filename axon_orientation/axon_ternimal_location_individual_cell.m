% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\phasemap_v1.1.1\phasemap'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'));
%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
folder = 'E:\BIL\download.brainimagelibrary.org+8811\0f\cd\0fcde5fdd6f7ccb2';
cd(folder)
dataDir = dir(folder);
for ii = 3:numel(dataDir)
    swcFileDir{ii-2} = dir(fullfile(folder, dataDir(ii).name, '*reg.swc')); % only look for files registered to Allen CCF
end
%% loading coords of each neuron from each mouse
ctr=0;
for mn = 1:numel(swcFileDir)
    for ii = 1:numel(swcFileDir{mn})
        ctr = ctr+1;
        thisFile = fullfile(swcFileDir{mn}(ii).folder, swcFileDir{mn}(ii).name);
        T = readtable(thisFile,'FileType','text');
        allCoords{ctr} = uint16(round([T.type T.x T.y T.z T.x__n T.parent])); % only keeping x y z coords and type of neuron 
        clear T thisFile
    end
end
%%
region_label = {"SSp-ll"};
st_region_indx = get_region_index_list(region_label,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
%%
% axon_region_label = {"Isocortex"};
axon_region_label = {"SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un","VISp"}; % axons in these areas
axon_st_cortex_indx = get_region_index_list(axon_region_label,st); % get all index values that contain "Isocortex" in st table
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%% plot all axons in a region, selected layers
h1 = figure;
clear region_all
regions_all1 = {"SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un"};
layer = {"1","2/3","4","5","6a","6b"};
regions_all = {};
for i = 1:6
% for i = kk
    ilayer = layer{i};
    region_test = cellfun(@(x) {append(x,ilayer)},regions_all1);
    regions_all = cat(2,regions_all,region_test);
end

cell_limit = 5;
[soma_all,axon_all] = get_all_axon_in_area_by_neuron(st,av,allCoords,regions_all,axon_st_cortex_indx);
%
clear vector_all S_diag angle1 polarity axon_current axon_current_2d soma_current_2d soma_center U S V
scale = 1;
current_axon_all = axon_all;
current_soma = soma_all;
cell_n = numel(current_axon_all);  
for icell = 1:cell_n
    axon_current = current_axon_all{icell};
    axon_current_2d = double(axon_current(:,[3,1]));
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
    soma_current_2d = current_soma{icell};
    soma_current_2d = soma_current_2d(1,:);
    soma_center(icell,:) = double(soma_current_2d(:,[3,1]));
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag(:,icell) = diag(S);
    vector_all(icell,:) = [U(1,1),U(2,1)];
    angle1(icell) = round(rad2deg(atan(U(2,1)/U(1,1))))+1;
    % plot(soma_center(1),soma_center(2),'.','MarkerFaceColor',color2(angle1,:),'MarkerSize',16)
end
% polarity = S_diag(1,:)./S_diag(2,:);
polarity = ones(size(vector_all,1),1);
%
color2 = hsv(180);
color3 = color2(angle1+90,:);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([100,500]);
ylim([300,800]);
scale1 = 10;
for i = 1:cell_n
    hold on;
    plot([soma_center(i,1)-vector_all(i,1)*scale1*polarity(i),soma_center(i,1)+vector_all(i,1)*scale1*polarity(i)],...
        [soma_center(i,2)-vector_all(i,2)*scale1*polarity(i),soma_center(i,2)+vector_all(i,2)*scale1*polarity(i)],...
        'color',color3(i,:),'LineWidth',1);
end
% scatter(soma_center(:,1),soma_center(:,2),[],color3,'filled');
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
title('All Layers');
%% plot all axons in a region by layer
% h1 = figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
h1 = figure;
for kk = 1:6
    clear region_all
    regions_all1 = {"SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un"};
    layer = {"1","2/3","4","5","6a","6b"};

    regions_all = {};
    % for i = 1:6
    for i = kk
        ilayer = layer{i};
        region_test = cellfun(@(x) {append(x,ilayer)},regions_all1);
        regions_all = cat(2,regions_all,region_test);
    end

    cell_limit = 5;
    % for i = 1:7
    %     % region_label = {"SSp-ll"};
    %     region_label = regions_all{i};
    %     % [soma_all{i},axon_all{i}] = get_all_axon_in_area(st,av,allCoords,region_label,axon_st_cortex_indx,cell_limit);
    %     % [soma_all{i},axon_all{i}] = get_all_axon_in_area(st,av,allCoords,region_label,axon_st_cortex_indx,cell_limit);
    %     [soma_all{i},axon_all{i}] = get_all_axon_in_area_by_neuron(st,av,allCoords,region_label,axon_st_cortex_indx);
    % end
    %
    [soma_all,axon_all] = get_all_axon_in_area_by_neuron(st,av,allCoords,regions_all,axon_st_cortex_indx);
    %
    clear vector_all S_diag angle1 polarity axon_current axon_current_2d soma_current_2d soma_center U S V
    scale = 1;
    current_axon_all = axon_all;
    current_soma = soma_all;
    cell_n = numel(current_axon_all);  
    for icell = 1:cell_n
        axon_current = current_axon_all{icell};
        axon_current_2d = double(axon_current(:,[3,1]));
        axon_center = mean(axon_current_2d,1);
        axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
        soma_current_2d = current_soma{icell};
        soma_current_2d = soma_current_2d(1,:);
        soma_center(icell,:) = double(soma_current_2d(:,[3,1]));
        [U,S,V] = svd(double(axon_current_2d)',"econ");
        S_diag(:,icell) = diag(S);
        vector_all(icell,:) = [U(1,1),U(2,1)];
        angle1(icell) = round(rad2deg(atan(U(2,1)/U(1,1))))+1;
        % plot(soma_center(1),soma_center(2),'.','MarkerFaceColor',color2(angle1,:),'MarkerSize',16)
    end
    polarity = S_diag(1,:)./S_diag(2,:);
    %
    color2 = hsv(180);
    % color2 = phasemap(180);
    color3 = color2(angle1+90,:);
    
    subplot(2,3,kk)
    overlayOutlines(coords,scale,[0.8,0.8,0.8]);
    set(gca, 'YDir','reverse');
    axis off; axis image;
    xlim([100,500]);
    ylim([400,800]);
    for i = 1:cell_n
        hold on;
        plot([soma_center(i,1)-vector_all(i,1)*10*polarity(i),soma_center(i,1)+vector_all(i,1)*10*polarity(i)],...
            [soma_center(i,2)-vector_all(i,2)*10*polarity(i),soma_center(i,2)+vector_all(i,2)*10*polarity(i)],...
            'color',color3(i,:),'LineWidth',1)
    end
    % scatter(soma_center(:,1),soma_center(:,2),[],color3,'filled');
    cb = colorbar;
    colormap(flipud(color2));
    cb.Ticks =  [0,0.5,1];
    cb.TickLabels = {'-90','0','90'};
    title(strcat("layer ", layer{kk}));
end
%% examples
scale = 1;
color1 = hsv(7);
for current_area = 1
    % current_axon_all = axon_all{current_area};
    % soma_current = soma_all{current_area};
    current_axon_all = axon_all;
    current_soma = soma_all;
    cell_n = numel(current_axon_all);  
    %%
    color2 = hsv(180);
    figure;
    % for icell = 1:cell_n
    batch = 80*0;
    for icell = 1:80
        axon_current = current_axon_all{icell+batch};
        axon_current_2d = double(axon_current(:,[3,1]));
        axon_center = mean(axon_current_2d,1);
        axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
        soma_current_2d = current_soma{icell+batch};
        soma_center = double(soma_current_2d(:,[3,1]));
        [U,S,V] = svd(double(axon_current_2d)',"econ");
        
        % subplot(1,cell_n,icell)
        subplot(10,8,icell)
        overlayOutlines(coords,scale,[0.8,0.8,0.8]);
        set(gca, 'YDir','reverse');
        hold on;
        plot(axon_current(:,3), axon_current(:,1), '.', 'Color', 'k','MarkerSize',2); 
        hold on;
        angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
        plot([soma_center(1)-U(1,1)*50,soma_center(1)+U(1,1)*50],[soma_center(2)-U(2,1)*50,soma_center(2)+U(2,1)*50],'color',color2(angle1,:),'LineWidth',2)
        axis off; axis image;
%         xlim([axon_center(1)-100,axon_center(1)+100]);
%         ylim([axon_center(2)-100,axon_center(2)+100]);
        xlim([100,500]);
        ylim([400,800]);
        % title(regions_all{current_area});
    end
end
%%
region_label = {"RSP"};
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
for i = 1:35
    ax1 = subplottight(7,5,i)
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