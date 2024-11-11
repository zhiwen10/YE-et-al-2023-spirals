function [hr2,p] = plotAxonOrientationMO4(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
load(fullfile(data_folder,'Revision','axons','MO_border_line_roi.mat'));
border_roi = roi;
load(fullfile(data_folder,'Revision','axons','MO_line_roi.mat'));
line_roi = roi;
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
center = [377,428]; % MOp
center2 = [244,542]; % SSp-un
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
%%
hr2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
ax2 = subplot(1,4,1);
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
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
scatter(center(1),center(2),36,'*','k');
hold on;
scatter(center2(1),center2(2),36,'*','k');
border_positions = border_roi.Position;
plot(border_positions(:,1),border_positions(:,2),'--','color',[0.5,0.5,0.5],'lineWidth',2);
% positions = roi.Position;
% plot(positions(:,1),positions(:,2),'k--','lineWidth',2);
xlim([0,600]);
ylim([0,1200]);

ax3 = subplot(1,4,2);
tf = inROI(roi,T1.soma_center_1,T1.soma_center_2);
T1 = T1(tf,:);
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
xlim([0,600]);
ylim([0,1200]);
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
scatter(center(1),center(2),36,'*','k');
% positions = roi.Position;
% positions(end+1,:) = positions(1,:);
% plot(positions(:,1),positions(:,2),'k','lineWidth',2);
positions = line_roi.Position;
plot(positions(:,1),positions(:,2),'k--','lineWidth',2);
%
ax4 =subplot(1,4,3);
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
center = [377,428]; % MOp
orthog_vector = soma_center1-center;
for k = 1:100
    rp_indx = randperm(size(axon_bias,1));
    vector_all1 = axon_bias(rp_indx,:);
    for i =1:size(axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
ax4 = plotAngleStairsShuffle(ax4,edges,N1,N_rp);
hold on;
xline(mean_bias,'r');
xlim([0,180]); 
ylim([5,25]);
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};
%%
ax3 = subplot(1,4,4);
Tq1 = T1(T1.soma_center_1<center(1) &T1.soma_center_2>center(2),:);
Tq2 = T1(T1.soma_center_1<center(1) &T1.soma_center_2<center(2),:);
Tq3 = T1(T1.soma_center_1>center(1) &T1.soma_center_2>center(2),:);
Tq4 = T1(T1.soma_center_1>center(1) &T1.soma_center_2<center(2),:);
colora3 = interp1(x1,color2,Tq1.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [Tq1.soma_center_1,Tq1.soma_center_2];
axon_bias = [Tq1.axon_bias_1,Tq1.axon_bias_2];
cell_n = size(Tq1,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*Tq1.pc_ratio(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*Tq1.pc_ratio(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*Tq1.pc_ratio(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*Tq1.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1),center(2),36,'*','k');
%%
[angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStatsMO(T1);
[h mu ul ll] = circ_mtest(angle_diff_real,deg2rad(90));
mean_rad = [mu,ll,ul];
mean_deg = rad2deg(mean_rad);
% for i = 1:100
for i = 1
    [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end
%%
print(hr2, fullfile(save_folder,'FigR2_circularbias_MO.pdf'),...
    '-dpdf', '-bestfit', '-painters');