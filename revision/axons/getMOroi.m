function getMOroi(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
% examples
h2ac = figure('Renderer', 'painters', 'Position', [100 100 400 400]);
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
ax1 = subplot(1,1,1);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
cell_n = size(T1,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
% xlim([0,600]);
% ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*T1.pc_ratio(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*T1.pc_ratio(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*T1.pc_ratio(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*T1.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(377,428,36,'*','k');
roi = drawpolygon(ax1);
save(fullfile(save_folder,'MO_roi.mat'),'roi');