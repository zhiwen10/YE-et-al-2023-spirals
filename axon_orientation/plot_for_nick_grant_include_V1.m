%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();
%%
% folder = 'E:\BIL\download.brainimagelibrary.org+8811\0f\cd\0fcde5fdd6f7ccb2';
% cd(folder)
% dataDir = dir(folder);
% for ii = 3:numel(dataDir)
%     swcFileDir{ii-2} = dir(fullfile(folder, dataDir(ii).name, '*reg.swc')); % only look for files registered to Allen CCF
% end
%% loading coords of each neuron from each mouse
% ctr=0;
% for mn = 1:numel(swcFileDir)
%     for ii = 1:numel(swcFileDir{mn})
%         ctr = ctr+1;
%         thisFile = fullfile(swcFileDir{mn}(ii).folder, swcFileDir{mn}(ii).name);
%         T = readtable(thisFile,'FileType','text');
%         allCoords{ctr} = uint16(round([T.type T.x T.y T.z T.x__n T.parent])); % only keeping x y z coords and type of neuron 
%         clear T thisFile
%     end
% end
%%
load('all_cell_with_parents.mat');
%%
region_label = {"SSp-ll"};
st_region_indx = get_region_index_list(region_label,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
%%
axon_region_label = {"VIS","RSP","AUD","VISa","VISC","VISrl","SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un","SSs"};
% axon_region_label = {};% axons in these areas
axon_st_cortex_indx = get_region_index_list(axon_region_label,st); % get all index values that contain "Isocortex" in st table

axon_region_label2 = {"Isocortex"};
axon_st_cortex_indx2 = get_region_index_list(axon_region_label2,st); % get all index values that contain "Isocortex" in st table
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
h1 = figure;
scale = 1;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis image; axis off;
print(h1, 'surface_outline', '-dpdf', '-bestfit', '-painters');
%% plot all axons in a region, selected layers
clear region_all
regions_all1 = {"VIS","RSP","AUD","VISa","VISC","VISrl","SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un","SSs"};
% regions_all1 = {"SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un"};
st_region_indx = get_region_index_list(regions_all1,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
type = 2; % only get axon
for i = 1:numel(cell_id)
    clear noi axonCoords axonCoords_filt
    axonCoords_filt = [];
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    axonCoords = noi(noi(:,1)== type,:);   
    % filter axon terminal in cortex
    % axon_st_cortex_indx = get_region_index_list(region_label,st); % get all index values that contain "Isocortex" in st table
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx2);
    axon_points_all{i} = axonCoords_filt;
end
[soma_all,axon_terminal_all] = get_all_axon_in_area_by_neuron(st,av,allCoords,regions_all1,axon_st_cortex_indx);
%%
h1 = figure;
clear vector_all S_diag angle1 polarity axon_current axon_current_2d soma_current_2d soma_center U S V
scale = 1;
a = cellfun(@(x) not(isempty(x)), axon_terminal_all);
axon_terminal_all = axon_terminal_all(a);
soma_all = soma_all(a);
cell_n = numel(axon_terminal_all);  
for icell = 1:cell_n
    axon_current = axon_terminal_all{icell};
    axon_current_2d = double(axon_current(:,[3,1]));
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
    soma_current_2d = soma_all{icell};
    soma_current_2d = soma_current_2d(1,:);
    soma_center(icell,:) = double(soma_current_2d(:,[3,1]));
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag(:,icell) = diag(S);
    vector_all(icell,:) = [U(1,1),U(2,1)];
    angle1(icell) = round(rad2deg(atan(U(2,1)/U(1,1))))+1;
    % plot(soma_center(1),soma_center(2),'.','MarkerFaceColor',color2(angle1,:),'MarkerSize',16)
end
polarity = S_diag(1,:)./S_diag(2,:)-1;
% polarity = ones(size(vector_all,1),1);
%
% color2 = hsv(180);
color2= colorcet('C06','N',180);
angle1(angle1==91) = 90;
color3 = color2(angle1+90,:);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
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

%% examples
scale = 1;
% color1 = hsv(7);
color1= colorcet('C06','N',7);
% current_axon_all = axon_all{current_area};
% soma_current = soma_all{current_area};
cell_n = numel(axon_terminal_all);  
%
% color2 = hsv(180);
color2= colorcet('C06','N',180);
h = figure;
% for icell = 1:cell_n
batch = 80*0;
% 61,58,43,59,98,138,193,50,122,51
% cells = [61,58,43,59];
cells = [155,145,86,131,98,207];
subplot(1,2,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');

subplot(1,2,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');

color4= colorcet('C06','N',numel(cells));

for icell = 1:numel(cells)    
    axon_current_all = axon_points_all{cells(icell)};
    axon_current = axon_terminal_all{cells(icell)};
    axon_current_2d = double(axon_current(:,[3,1]));
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);
    soma_current_2d = soma_all{cells(icell)};
    soma_center1 = double(soma_current_2d(:,[3,1]));
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;

    subplot(1,2,1);
%     overlayOutlines(coords,scale,[0.8,0.8,0.8]);
%     set(gca, 'YDir','reverse');
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color2(angle1,:),'MarkerSize',2); 
    hold on;
    scatter(soma_center1(1),soma_center1(2),10,'k','filled');
    axis off; axis image;
    xlim([100,500]);
    ylim([350,800]);
    
    subplot(1,2,2);
%     overlayOutlines(coords,scale,[0.8,0.8,0.8]);
%     set(gca, 'YDir','reverse');
    hold on;
    scatter(axon_current(:,3), axon_current(:,1),6,'MarkerFaceColor',...
        color2(angle1,:),'MarkerEdgeColor','None','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1); 
    hold on;
    plot([soma_center1(1)-U(1,1)*30*polarity3(icell),soma_center1(1)+U(1,1)*30*polarity3(icell)],...
        [soma_center1(2)-U(2,1)*30*polarity3(icell),soma_center1(2)+U(2,1)*30*polarity3(icell)],'color',color2(angle1,:),'LineWidth',2)
    hold on;
    scatter(soma_center1(1),soma_center1(2),10,'k','filled');
    axis off; axis image;
    hold on;
    xlim([100,500]);
    ylim([350,800]);
end


%% examples
scale = 1;
color1= colorcet('C06','N',7);
cell_n = numel(axon_terminal_all);  
color2= colorcet('C06','N',180);
h = figure;
subplot(1,1,1)
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for i = 1:cell_n
    hold on;
    plot([soma_center(i,1)-vector_all(i,1)*30*polarity(i),soma_center(i,1)+vector_all(i,1)*30*polarity(i)],...
        [soma_center(i,2)-vector_all(i,2)*30*polarity(i),soma_center(i,2)+vector_all(i,2)*30*polarity(i)],...
        'color',color3(i,:),'LineWidth',1)
end
axis off; axis image;
xlim([0,600]);
ylim([300,1100]);
set(gca,'Visible',false)
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
%%
print(h, 'cortex_axon_lateral_bias3_cetc06', '-dpdf', '-bestfit', '-painters');
%%
center = [244,542]; % SSp-un
orthog_vector = soma_center-center;
%%
for i =1:size(vector_all,1)
    u = [vector_all(i,:) 0];
    v = [orthog_vector(i,:) 0];
    ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
end    
ang_diff2 = rad2deg(ang_diff1);
edges = 0:10:180;
[N1,edges] = histcounts(ang_diff2,edges);
%% shuffle
for k = 1:100
    rp_indx = randperm(size(vector_all,1));
    vector_all1 = vector_all(rp_indx,:);
    for i =1:size(vector_all1,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
%%
h1 = figure;
stairs(edges(1:end-1),N1,'r');
hold on;
mean_N_rp = mean(N_rp,1);
std_N_rp = std(N_rp,[],1)/sqrt(size(N_rp,1));
stairs(edges(1:end-1),mean_N_rp+std_N_rp,'k');
hold on;
stairs(edges(1:end-1),mean_N_rp-std_N_rp,'k');
hold on;
stairs(edges(1:end-1),mean_N_rp,'b');

%Calculate the error bars
uE = mean_N_rp+std_N_rp;
lE = mean_N_rp-std_N_rp;
uEs = uE(1:end-1);
lEs = lE(1:end-1);
uE_new(1:2:2*numel(uE)-1) = uE;
uE_new(2:2:2*numel(uE)-1) = uEs;
lE_new(1:2:2*numel(lE)-1) = lE;
lE_new(2:2:2*numel(lE)-1) = lEs;

xu = edges(1:end-1);
xus = xu(2:end);
xu_new(1:2:2*numel(xu)-1) = xu;
xu_new(2:2:2*numel(xu)-1) = xus;

%Make the patch (the shaded error bar)
yP=[lE_new,fliplr(uE_new)];
xP=[xu_new,fliplr(xu_new)];
Hpatch=patch(xP,yP,1);
patchColor=[0.5,0.5,0.5]; faceAlpha = 0.5;
set(Hpatch,'facecolor',patchColor, ...
    'edgecolor','none', ...
    'facealpha',faceAlpha, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')

xticks([0 30 60 90 120 150 180])
xticklabels({'0' '30' '60' '90' '120' '150' '180'})

xlabel('angle bins');
ylabel('bin counts');
legend({'data','shuffle+std','shuffle-std','100X shuffle'});
%%
print(h1, 'axon_1pc_orthog_angle-1', '-dpdf', '-bestfit', '-painters')
%% examples
axon_points_all2 = {};
soma_all2 = {};
for k = 1:numel(regions_all1)
% k = 1;
    iregion = regions_all1{k};
    st_region_indx = get_region_index_list(iregion,st); % get all index values that contain "SSp-ll" in st table
    [cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
    [soma_all_temp,axon_terminal_all1] = get_all_axon_in_area_by_neuron(st,av,allCoords,iregion,axon_st_cortex_indx);
    type = 2; % only get axon
    type1 = 1;
    type2 = 4; % only get dendrite
    axon_points_all1 = {};
    soma_all1 = {};
    for i = 1:numel(cell_id)
        clear noi axonCoords axonCoords_filt
        axonCoords_filt = [];
        noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
        axonCoords1 = noi(noi(:,1)== type,:);   
        % filter axon terminal in cortex
        % axon_st_cortex_indx = get_region_index_list(region_label,st); % get all index values that contain "Isocortex" in st table
        hemi =2;
        [axonCoords_filt] = filter_axon_position(axonCoords1,av,axon_st_cortex_indx2,hemi);
        axon_points_all1{i} = axonCoords_filt;
        %%
        soma = noi(noi(:,1)== type1,2:4)/10;  
        soma_all1{i} = soma;
        %%
        dendrite = noi(noi(:,1)== type2,2:4)/10;  
        dendrite_all1{i} = dendrite;
    end
    axon_points_all2{k} = axon_points_all1;
    soma_all2{k} = soma_all1;
    dendrite_all2{k} = dendrite_all1;
    axon_terminal_all2{k} = axon_terminal_all1;
end

%%
cell_n = cellfun(@(x) numel(x),axon_points_all2);
%%
regions_all1 = {"VIS","RSP","AUD","VISa","VISC","VISrl","SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un","SSs"};
color4 = hsv(7);
for k=2
    current_axon = axon_points_all2{k};
    current_soma = soma_all2{k};
    current_dendrite = dendrite_all2{k};
    current_terminal = axon_terminal_all2{k};
    cell_n = numel(current_axon);
    if cell_n>40
        cell_n =40;
    end
    h = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
    cell_id = 1:cell_n;
    for icell = 1:numel(cell_id)    
        clear axon_current_all
        axon_current_all = current_axon{cell_id(icell)};
        current_soma_all = current_soma{cell_id(icell)};
        current_dendrite_all = current_dendrite{cell_id(icell)};
        subplot(8,5,icell);
        overlayOutlines(coords,scale,[0.8,0.8,0.8]);
        set(gca, 'YDir','reverse');
        hold on;
        scatter(axon_current_all(:,3), axon_current_all(:,1),0.5,'MarkerFaceColor',color4(k,:),'MarkerEdgeColor','None'); 
        % plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color4(k,:),'MarkerSize',0.5); 
        hold on;
        scatter(current_soma_all(:,3),current_soma_all(:,1),2,'k','filled');
        hold on;
        scatter(current_dendrite_all(:,3),current_dendrite_all(:,1),2,'r','filled');
        axis off; axis image;
    %     xlim([100,500]);
    %     ylim([350,800]);
    end
    savefig(h,[char(regions_all1{k}) '-' num2str(cell_id(1))]);
    print(h, [char(regions_all1{k}) '-' num2str(cell_id(1))], '-dpdf', '-bestfit', '-painters');
end
%%
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
% 3 subareas
root1 = {'/997/'};
maskPath{1}  = '/997/8/343/1129/549/'; % Th
maskPath{2}  = '/997/8/567/623/477/'; % STR
maskPath{3}  = '/997/8/343/313/'; % midbrain
maskPath{4}  = "/997/8/567/688/695/315/"; % isocortex
maskPath{5}  = "/997/8/567/688/695/698/"; % olfactory

atlas1{1} = double(squeeze(atlas(80,:,:)))';
atlas1{2} = double(squeeze(atlas(:,:,190)));
atlas1{3} = double(squeeze(atlas(:,160,:)));

% order_2d(:,1) = [3,1];
% order_2d(:,2) = [3,2];
% order_2d(:,3) = [1,2];

order_2d(:,1) = [1,3];
order_2d(:,2) = [1,2];
order_2d(:,3) = [3,2];
axis_lim(:,1) = [264,264];
axis_lim(:,2) = [264,264];
axis_lim(:,3) = [264,264];
%%
regions_all1 = {"VIS","RSP","AUD","VISa","VISC","VISrl","SSp-tr","SSp-ll","SSp-ul","SSp-m","SSp-n","SSp-bfd","SSp-un","SSs"};
color4 = hsv(14);
for k=2
    current_axon = axon_points_all2{k};
    current_soma = soma_all2{k};
    current_dendrite = dendrite_all2{k};
     current_terminal = axon_terminal_all2{k};
    cell_n = numel(current_axon);
    if cell_n>40
        cell_n =40;
    end
    for kk = 1:2
        h = figure('Renderer', 'painters', 'Position', [100 100 700 900]);
        cell_id = 1:cell_n;
        count1 = 1;
        % for icell = 1:numel(cell_id)    
        for icell = (kk-1)*8+1:(kk-1)*8+8
            if icell<cell_n
            clear axon_current_all
            axon_current_all = current_axon{cell_id(icell)}/5;
            current_soma_all = current_soma{cell_id(icell)}/5;
            current_dendrite_all = current_dendrite{cell_id(icell)}/5;
            for iplot = 1:3
                subplot(8,3,(count1-1)*3+iplot);
                scale3 = 1;
                plotOutline_fill(root1,st,atlas1{iplot},[],scale3,'w');
    %             plotOutline_fill(maskPath(1),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
    %             plotOutline_fill(maskPath(2),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
    %             plotOutline_fill(maskPath(3),st,atlas1{iplot},[],scale3,[0.5,0.5,0.5]);
                if (iplot == 2) | (iplot==3)
                    plotOutline_fill(maskPath(4:5),st,atlas1{iplot},[],scale3,'w');
                end
                axis image; axis off;
                set(gca, 'YDir','reverse');
                xlim([-20,axis_lim(1,iplot)]+20);
                ylim([-20,axis_lim(2,iplot)]+20);
                hold on;
                scatter(axon_current_all(:,order_2d(1,iplot)), axon_current_all(:,order_2d(2,iplot)),0.5,'MarkerFaceColor',color4(k,:),'MarkerEdgeColor','None'); 
                % plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color4(k,:),'MarkerSize',0.5); 
                hold on;
                scatter(current_soma_all(:,order_2d(1,iplot)),current_soma_all(:,order_2d(2,iplot)),2,'k','filled');
    %             hold on;
    %             scatter(current_dendrite_all(:,order_2d(1,iplot)),current_dendrite_all(:,order_2d(2,iplot)),2,'r','filled');
                axis off; axis image;
            end
            count1 = count1+1;
        %     xlim([100,500]);
        %     ylim([350,800]);
            end
        end
        saveas(h,[char(regions_all1{k}) '-' num2str(cell_id(1)) '-' num2str(kk) '.png']);
        print(h, [char(regions_all1{k}) '-' num2str(cell_id(1)) '-' num2str(kk)], '-dpdf', '-bestfit', '-painters');
    end
end
%%
axon_current_all1{7}  = current_axon{2};
soma_current_2d1{7} = current_soma{2};
axon_current1{7} = current_terminal{2};
%%
axon_current_all1{8}  = current_axon{15};
soma_current_2d1{8} = current_soma{15};
axon_current1{8} = current_terminal{15};
%% examples
scale = 1;
color1= colorcet('C06','N',8);
cell_n = numel(axon_terminal_all);  
color2= colorcet('C06','N',180);
h = figure;
subplot(1,4,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');

subplot(1,4,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
color4= colorcet('C06','N',numel(cells));

for icell = 1:8
    axon_current_all = axon_current_all1{icell};
    axon_current = axon_current1{icell};
    soma_current_2d = soma_current_2d1{icell};
    
    axon_current_2d = double(axon_current(:,[3,1]));
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);    
    soma_center1 = double(soma_current_2d(:,[3,1]));
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;

    subplot(1,4,1);
%     overlayOutlines(coords,scale,[0.8,0.8,0.8]);
%     set(gca, 'YDir','reverse');
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color2(angle1,:),'MarkerSize',2); 
    hold on;
    scatter(soma_center1(1),soma_center1(2),10,'k','filled');
    axis off; axis image;
    xlim([0,600]);
    ylim([0,1200]);
    
    subplot(1,4,2);
%     overlayOutlines(coords,scale,[0.8,0.8,0.8]);
%     set(gca, 'YDir','reverse');
    hold on;
    scatter(axon_current(:,3), axon_current(:,1),6,'MarkerFaceColor',...
        color2(angle1,:),'MarkerEdgeColor','None','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1); 
    hold on;
    plot([soma_center1(1)-U(1,1)*30*polarity3(icell),soma_center1(1)+U(1,1)*30*polarity3(icell)],...
        [soma_center1(2)-U(2,1)*30*polarity3(icell),soma_center1(2)+U(2,1)*30*polarity3(icell)],'color',color2(angle1,:),'LineWidth',2)
    hold on;
    scatter(soma_center1(1),soma_center1(2),10,'k','filled');
    axis off; axis image;
    hold on;
    xlim([0,600]);
    ylim([0,1200]);
end

color1= colorcet('C06','N',7);
cell_n = numel(axon_terminal_all);  
color2= colorcet('C06','N',180);
subplot(1,4,3)
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for i = 1:cell_n
    hold on;
    plot([soma_center(i,1)-vector_all(i,1)*30*polarity(i),soma_center(i,1)+vector_all(i,1)*30*polarity(i)],...
        [soma_center(i,2)-vector_all(i,2)*30*polarity(i),soma_center(i,2)+vector_all(i,2)*30*polarity(i)],...
        'color',color3(i,:),'LineWidth',1)
end
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
subplot(1,4,3);
set(gca,'Visible',false)
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-90','0','90'};
axesPosition = get(cb, 'Position');
cb.Position([1,2,3,4]) = [0.8,0.3,0.01,0.4];
%%
print(h, 'all_cortex_axon_lateral_bias3_cetc06', '-dpdf', '-bestfit', '-painters');