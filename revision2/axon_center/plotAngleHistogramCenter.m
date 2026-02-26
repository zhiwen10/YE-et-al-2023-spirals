function h2d = plotAngleHistogramCenter(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
% T = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells_ssp.csv'));
T = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells.csv'));
soma_center = [T.soma_center_1,T.soma_center_2];
axon_vector = [T.axon_bias_1,T.axon_bias_2];
edges = 0:10:180;
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
[centers_x,centers_y] = meshgrid(0:10:1140,0:10:1320);
centers_x1 = centers_x(:);
centers_y1 = centers_y(:);
clear N1;
centerN = numel(centers_x1);
for k = 1:centerN
    ang_diff1 = zeros(size(axon_vector,1),1);
    center = [centers_x1(k),centers_y1(k);];
    orthog_vector = soma_center-center;
    for i =1:size(axon_vector,1)
        u = [axon_vector(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i,1) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N1(k,:),edges] = histcounts(ang_diff2,edges);
end
N1 = reshape(N1,size(centers_x,1),size(centers_x,2),18);
%%
% centers2(1,:) = [220,400];
% centers2(2,:) = [220,560];
% centers2(3,:) = [220,800];
centers2(1,:) = [260,480];
centers2(2,:) = [260,560];
centers2(3,:) = [270,670];
for i = 1:3
    indx(i,1) = find(centers_x1 == centers2(i,1) & centers_y1 == centers2(i,2));
end
%%
N2 =  reshape(N1, size(N1,1)*size(N1,2),size(N1,3));
% peaks = squeeze(mean(N1(:,:,8:10),3));
peaks = squeeze(mean(N1(:,:,6:12),3));
troughs = squeeze(mean(N1(:,:,[1,2,3,16,17,18]),3));
% ratio = peaks./troughs;
ratio = peaks-troughs;
ratio1 = imresize(ratio,[1320,1140]);
Ni =ratio(indx);
%%
[a1,b1] = ismember([centers_x1,centers_y1],brain_index,'rows');
%%
color1 = colormap(hot(100));
mina = min(ratio(:));
maxa = max(ratio(:));
scale1 = linspace(mina,maxa,100);
color3 = interp1(scale1,color1,ratio(:));
centers_x2 = centers_x1;
centers_y2 = centers_y1;
centers_x2(not(a1)) = nan;
centers_y2(not(a1)) = nan;
color4 = color3;
color4(not(a1),:) = nan;
%%
BW2 = logical(projectedAtlas1);
h2 = figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);

ax1 = subplot(3,3,[1,4,7]);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T.soma_center_1,T.soma_center_2];
axon_bias = [T.axon_bias_1,T.axon_bias_2];
cell_n = size(T,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
% axis off; 
axis image;
xlim([0,600]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*T.pc_ratio(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*T.pc_ratio(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*T.pc_ratio(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*T.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(centers2(:,1), centers2(:,2),36,'*','k');
cb = colorbar;
colormap(ax1,flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-pi/2','0','pi/2'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.04 a(2)+0.05 a(3) a(4)/2])% To change size

color1 = {'r','g','k'};
for i = 1:3    
    ax1 = subplot(3,3,2+3*(i-1));
    stairs(ax1,edges(1:end-1),squeeze(N2(indx(i),:)),'color',color1{i});
    xlim([0,180]); 
    ylim([0,40]);
    xticks([0 45 90 135 180])
    xticklabels({'0','1/4*pi','1/2*pi','3/4*pi','pi'});   
    box off;
    text(90,35,['Center index = ' num2str(round(Ni(i)))]);
end

ax2 = subplot(3,3,[3,6,9]);
im1 = imagesc(ones(size(ratio1)));
hold on;
scatter(centers_x2,centers_y2,3,color4,'filled');
hold on;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
% axis off; 
axis image;
set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
scatter(centers2(:,1), centers2(:,2),36,'*','k');
colormap(ax2, hot);
caxis([0,1]);
cb2 = colorbar;
cb2.Ticks = linspace(0,1,4);
cb2.TickLabels = {'-10','0','10','20'};
a =  cb2.Position; %gets the positon and size of the color bar
set(cb2,'Position',[a(1)+0.04 a(2)+0.05 a(3) a(4)/2])% To change size
xlim([0,600]);
ylim([0,1200]);
%%
print(h2,'axon_center3.pdf','-dpdf', '-bestfit', '-painters');