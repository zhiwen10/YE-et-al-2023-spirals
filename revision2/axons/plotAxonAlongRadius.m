function [h2ac] = plotAxonAlongRadius(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells.csv'));
load(fullfile(data_folder, 'axons','example_12_cells.mat'));
%% examples
h2ac = figure('Renderer', 'painters', 'Position', [100 100 600 700]);
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
% center = [262,597]; % SSp-un
% center = [236,629]; % SSp-un
% center = [278,616]; % SSp-un
center = [244,542]; % SSp-un

ax1 = subplot(2,3,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
scatter(center(1), center(2),36,'*','k');
hold on;
for icell = [1:8,10:12]
    axon_current_all = axon_current_all1{icell};
    axon_current = axon_current1{icell};
    soma_current_2d = soma_current_2d1{icell};
    
    axon_current_2d = double(axon_current(:,[3,1]));
    
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-...
        repmat(axon_center,size(axon_current_2d,1),1);    
    soma_center1 = double(soma_current_2d(:,[3,1]));
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1),...
        '.', 'Color', color2(angle1,:),'MarkerSize',1); 
    hold on;
    scatter(axon_current(:,3),axon_current(:,1),1,[0.5,0.5,0.5],'filled');
    hold on;
    plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),...
        soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
        [soma_center1(2)-U(2,1)*scale1*polarity3(icell),...
        soma_center1(2)+U(2,1)*scale1*polarity3(icell)],...
        'color','k','LineWidth',1)
    axis off; axis image;
    hold on;
    xlim([0,600]);
    ylim([0,1200]);
end
hold on;
plot([200,400],[1100,1100],'k');
text(300,1150,'2 mm');

ax2 = subplot(2,3,2);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1.bias_angle);
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
scatter(center(1),center(2),'*','k');
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-pi/2','0','pi/2'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.04 a(2)+0.05 a(3) a(4)/2])% To change size

ax3 = subplot(2,3,3);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
cell_n = size(T1,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
for r = 50:50:200
    radian_bins = [0:pi/20:2*pi];
    x1 = center(1)+r*cos(radian_bins);
    y1 = center(2)+r*sin(radian_bins);
    % hold on;
    % scatter(x1,y1,2,'k');
    plot(x1,y1,'k');
    hold on;
end
for r = 200:50:350
    radian_bins = [-1/2*pi:pi/20:5/6*pi];
    x1 = center(1)+r*cos(radian_bins);
    y1 = center(2)+r*sin(radian_bins);
    % hold on;
    % scatter(x1,y1,2,'k');
    plot(x1,y1,'k');
    hold on;
end
%
ax4 = subplot(2,3,4);
% scatter(distance,T1.pc_ratio);
lm = fitlm(distance,T1.pc_ratio);
% hold on;
plot(lm);



ax5 = subplot(2,3,5);
coef1 = corrcoef(distance,T1.pc_ratio);
coef1 = coef1(1,2);
data_n = length(distance);
for i = 1:5000
    pa = randperm(data_n);
    distance1 = distance(pa);
    corr_temp = corrcoef(distance1,T1.pc_ratio);
    coef2(i,1) = corr_temp(1,2);
end
histogram(coef2);
hold on;
xline(coef1,'r--');
[h,p] = ttest2(coef2,coef1);

ax6 = subplot(2,3,6);
soma_ref = soma_center1-center;
distance = vecnorm(soma_ref,2,2);
radius_bins = 0:25:350;
pc_ratio_bins = {};
for i = 1:length(radius_bins)-1
    pc_ratio_index = (distance>radius_bins(i) & distance<radius_bins(i+1));
    pc_ratio_temp = T1.pc_ratio(pc_ratio_index);
    pc_ratio_bins{i} = pc_ratio_temp;
end
pc_ratio_bins_mean = cellfun(@mean, pc_ratio_bins);
pc_ratio_bins_std = cellfun(@std, pc_ratio_bins);
pc_ratio_bins_count = cellfun(@length,pc_ratio_bins);
pc_ratio_bins_sem = pc_ratio_bins_std./sqrt(pc_ratio_bins_count);
radius_bins2 = radius_bins(2:end);
% for i = 1:length(radius_bins)-1
%     pc_ratio_bins1 = pc_ratio_bins{i};
%     radius_temp = ones(size(pc_ratio_bins1))*radius_bins(i+1);
%     scatter(radius_temp,pc_ratio_bins1,6,'r');
%     hold on;
% end
index = (pc_ratio_bins_count>4);
errorbar(radius_bins2(index),pc_ratio_bins_mean(index), pc_ratio_bins_sem(index));
print(h2ac, fullfile(save_folder,'Fig2ac_axon_bias_all2.pdf'),...
    '-dpdf', '-bestfit', '-painters');