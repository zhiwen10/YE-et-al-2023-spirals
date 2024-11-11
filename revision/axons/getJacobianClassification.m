function [hr2] = plotAxonOrientationMO4(data_folder,save_folder)
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
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
tf = inROI(roi,T1.soma_center_1,T1.soma_center_2);
T1 = T1(tf,:);
Tq1 = T1(T1.soma_center_1<center(1) &T1.soma_center_2<center(2),:); % top left
Tq2 = T1(T1.soma_center_1>center(1) &T1.soma_center_2<center(2),:); % top right
Tq3 = T1(T1.soma_center_1<center(1) &T1.soma_center_2>center(2),:); %bottom left
Tq4 = T1(T1.soma_center_1>center(1) &T1.soma_center_2>center(2),:); % bottom right
%%
hr2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
ax2 = subplot(1,4,1);
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
positions = line_roi.Position;
plot(positions(:,1),positions(:,2),'k--','lineWidth',2);
%
ax3 = subplot(1,4,2);
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
for i = 1:1000
    indx1 = randperm(size(Tq1,1),1);
    indx2 = randperm(size(Tq2,1),1);
    indx3 = randperm(size(Tq3,1),1);
    indx4 = randperm(size(Tq4,1),1);
    vector1 = [Tq1.axon_bias_1(indx1),Tq1.axon_bias_2(indx1)];
    vector2 = [Tq2.axon_bias_1(indx2),Tq2.axon_bias_2(indx2)];
    vector3 = [Tq3.axon_bias_1(indx3),Tq3.axon_bias_2(indx3)];
    vector4 = [Tq4.axon_bias_1(indx4),Tq4.axon_bias_2(indx4)];
    vx = [vector1(1),vector2(1);vector3(1),vector4(1)];
    vy = [vector1(2),vector2(2);vector3(2),vector4(2)];
    itype{i,1} = getJacobian(vx,vy);
end    
itype = string(itype);
[C1, ia1, ic1] = unique(itype);
a_counts = accumarray(ic1,1);
value_counts = [C1, a_counts];
ratio = a_counts(C1 == 'nodes')./sum(a_counts);
%%
for k = 1:1000
    indxa = randperm(size(T1,1));
    T1a = T1;
    T1a.axon_bias_1 = T1a.axon_bias_1(indxa);
    T1a.axon_bias_2 = T1a.axon_bias_2(indxa);
    T1a.bias_angle = T1a.bias_angle(indxa);
    Tq1a = T1a(T1a.soma_center_1<center(1) &T1a.soma_center_2<center(2),:); % top left
    Tq2a = T1a(T1a.soma_center_1>center(1) &T1a.soma_center_2<center(2),:); % top right
    Tq3a = T1a(T1a.soma_center_1<center(1) &T1a.soma_center_2>center(2),:); %bottom left
    Tq4a = T1a(T1a.soma_center_1>center(1) &T1a.soma_center_2>center(2),:); % bottom right
    for i = 1:1000
        indx1 = randperm(size(Tq1a,1),1);
        indx2 = randperm(size(Tq2a,1),1);
        indx3 = randperm(size(Tq3a,1),1);
        indx4 = randperm(size(Tq4a,1),1);
        vector1 = [Tq1a.axon_bias_1(indx1),Tq1a.axon_bias_2(indx1)];
        vector2 = [Tq2a.axon_bias_1(indx2),Tq2a.axon_bias_2(indx2)];
        vector3 = [Tq3a.axon_bias_1(indx3),Tq3a.axon_bias_2(indx3)];
        vector4 = [Tq4a.axon_bias_1(indx4),Tq4a.axon_bias_2(indx4)];
        vx = [vector1(1),vector2(1);vector3(1),vector4(1)];
        vy = [vector1(2),vector2(2);vector3(2),vector4(2)];
        itype{i,1} = getJacobian(vx,vy);
    end    
    itype = string(itype);
    [C1, ia1, ic1] = unique(itype);
    a_counts = accumarray(ic1,1);
    value_counts = [C1, a_counts];
    ratio(k,1) = a_counts(C1 == 'nodes')./sum(a_counts);
end
%%
hr2 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
ax2 = subplot(1,4,1);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1a.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T1a.soma_center_1,T1a.soma_center_2];
axon_bias = [T1a.axon_bias_1,T1a.axon_bias_2];
cell_n = size(T1a,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*T1a.pc_ratio(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*T1a.pc_ratio(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*T1a.pc_ratio(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*T1a.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1),center(2),36,'*','k');
positions = line_roi.Position;
plot(positions(:,1),positions(:,2),'k--','lineWidth',2);
%%
print(hr2, fullfile(save_folder,'FigR2_circularbias_MO.pdf'),...
    '-dpdf', '-bestfit', '-painters');