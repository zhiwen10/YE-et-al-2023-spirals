function [hr3,p] = plotAxonOrientationMO3(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
center = [377,428]; % MOp
center2 = [244,542]; % SSp-un
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
%%
hr3 = figure('Renderer', 'painters', 'Position', [100 100 800 500]);
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO2.csv'));
ax1 = subplot(2,4,1);
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
xlim([0,600]);
ylim([0,1200]);

ax2 = subplot(2,4,2);
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
% positions = roi.Position;
% positions(end+1,:) = positions(1,:);
% plot(positions(:,1),positions(:,2),'k','lineWidth',2);
xlim([0,600]);
ylim([0,1200]);

ax3 = subplot(2,4,3);
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
positions = roi.Position;
positions(end+1,:) = positions(1,:);
plot(positions(:,1),positions(:,2),'k','lineWidth',2);
hold on;
plot([center(1),center(1)],[center(2)+200,center(2)-200],'k--','lineWidth',2);
hold on;
plot([center(1),center(1)-200],[center(2),center(2)],'k--','lineWidth',2);
%
ax4 =subplot(2,4,4);
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
center = [377,428]; % MOp
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
T3 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
tf = inROI(roi,T3.soma_center_1,T3.soma_center_2);
T3 = T3(tf,:);
axon_vector = [T3.axon_bias_1,T3.axon_bias_2];
theta = 90;
for i = 1:size(axon_vector,1)
    v1 = axon_vector(i,:);
    vr = v1*[1;1i]*exp(-1i*theta*pi/180);
    angle1(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
end
angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
angle2a = angle1(T3.soma_center_1<center(1) & T3.soma_center_2>center(2));
angle2b = angle1(T3.soma_center_1>center(1));
angle2c = angle1(T3.soma_center_1<center(1) & T3.soma_center_2<center(2));
%
ax5 =subplot(2,4,5);
hist2a = histogram(angle2a,10,'faceColor','None','EdgeColor','k','DisplayStyle','stairs');
angle2a_max = mean(hist2a.BinEdges(hist2a.Values == max(hist2a.Values)));
xline(angle2a_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
ax5.XTick =  [-90,-45,0,45,90];
ax5.XTickLabel = {'-pi/2','-pi/4','0','pi/4','pi/2'};

ax6 =subplot(2,4,6);
hold on;
hist2b = histogram(angle2b,10,'faceColor','None','EdgeColor','r','DisplayStyle','stairs');
angle2b_max = hist2b.BinEdges(hist2b.Values == max(hist2b.Values));
xline(angle2b_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
ax6.XTick =  [-90,-45,0,45,90];
ax6.XTickLabel = {'-pi/2','-pi/4','0','pi/4','pi/2'};

ax7 =subplot(2,4,7);
hold on;
hist2c = histogram(angle2c,10,'faceColor','None','EdgeColor','b','DisplayStyle','stairs');
angle2c_max = hist2c.BinEdges(hist2c.Values == max(hist2c.Values));
xline(angle2c_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
ax7.XTick =  [-90,-45,0,45,90];
ax7.XTickLabel = {'-pi/2','-pi/4','0','pi/4','pi/2'};
%%
angle_all = [angle2a_max+5,angle2b_max+5,angle2c_max+5];
angle_all2 = angle_all+90;
angle_all2(2:3) = angle_all2(2:3)-180;
color2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,angle_all2);
%
angle_all3 = angle_all;
angle_all3(2) = angle_all3(2)-90;
r = 1;
z = r*exp(1i*deg2rad(angle_all3));
z=z';
%
index2a = (T3.soma_center_1<center(1) & T3.soma_center_2>center(2));
soma2a = mean([T3.soma_center_1(index2a),T3.soma_center_2(index2a)],1);
index2b = (T3.soma_center_1>center(1));
soma2b = mean([T3.soma_center_1(index2b),T3.soma_center_2(index2b)],1);
index2c = (T3.soma_center_1<center(1) & T3.soma_center_2<center(2));
soma2c = mean([T3.soma_center_1(index2c),T3.soma_center_2(index2c)],1);
%%
hist2a.EdgeColor = colora3(1,:);
hist2b.EdgeColor = colora3(2,:);
hist2c.EdgeColor = colora3(3,:);
%%
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
ax8 = subplot(2,4,8);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [soma2a;soma2b;soma2c];
axon_bias = [real(z),imag(z)];
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
scale1 = 50;
for i = 1:3
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1,...
        soma_center1(i,1)+axon_bias(i,1)*scale1],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1,...
        soma_center1(i,2)+axon_bias(i,2)*scale1],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1),center(2),36,'*','k');
xlim([0,600]);
ylim([0,1200]);
positions = roi.Position;
positions(end+1,:) = positions(1,:);
plot(positions(:,1),positions(:,2),'k','lineWidth',2);
hold on;
plot([center(1),center(1)],[center(2)+200,center(2)-200],'k--','lineWidth',2);
hold on;
plot([center(1),center(1)-200],[center(2),center(2)],'k--','lineWidth',2);
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
print(hr3, fullfile(save_folder,'FigR3_circularbias_MO.pdf'),...
    '-dpdf', '-bestfit', '-painters');