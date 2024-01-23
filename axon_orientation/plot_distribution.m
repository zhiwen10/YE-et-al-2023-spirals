function h = plot_distribution(T1,coords)
%% examples
h = figure('Renderer', 'painters', 'Position', [100 100 1100 300]);
load('example_12_cells.mat');
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);

subplot(1,5,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
hold on;
center = [244,542]; % SSp-un
scatter(244,542,36,'*','k');
hold on;
plot([200,400],[1100,1100],'k');
text(300,1150,'2 mm');
v = [100,400;100,800;500,800;500,400]/scale;
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
axis off; 
axis image;
xlim([0,600]);
ylim([0,1200]);

subplot(1,5,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for icell = 1
    axon_current_all = axon_current_all1{icell};
    axon_current = axon_current1{icell};
    soma_current_2d = soma_current_2d1{icell};
    
    axon_current_2d = double(axon_current(:,[3,1]));
    
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-repmat(axon_center,size(axon_current_2d,1),1);    
    soma_center1 = double(soma_current_2d(:,[3,1]));
    sample_n = size(axon_current_2d,1);
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color2(angle1,:),'MarkerSize',1); 
    hold on;
    scatter(axon_current(:,3),axon_current(:,1),1,[0.5,0.5,0.5],'filled');
    hold on;
%     plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
%         [soma_center1(2)-U(2,1)*scale1*polarity3(icell),soma_center1(2)+U(2,1)*scale1*polarity3(icell)],'color','k','LineWidth',1);
    plot([soma_center1(1)-U(1,1)*S_diag3(1)/sqrt(sample_n),soma_center1(1)+U(1,1)*S_diag3(1)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,1)*S_diag3(1)/sqrt(sample_n),soma_center1(2)+U(2,1)*S_diag3(1)/sqrt(sample_n)],'color','k','LineWidth',1);
    hold on;
    plot([soma_center1(1)-U(1,2)*S_diag3(2)/sqrt(sample_n),soma_center1(1)+U(1,2)*S_diag3(2)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,2)*S_diag3(2)/sqrt(sample_n),soma_center1(2)+U(2,2)*S_diag3(2)/sqrt(sample_n)],'color','k','LineWidth',1);
    axis off; 
    axis image;
end
hold on;
center = [244,542]; % SSp-un
scatter(244,542,36,'*','k');
hold on;
plot([soma_center1(1),244],[soma_center1(2),542],'--','color','k');
hold on;
plot([150,250],[800,800],'k');
text(150,850,'1 mm');
xlim([100,500]);
ylim([400,800]);

subplot(1,5,3);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
scatter(244,542,36,'*','k');
hold on;
for icell = [1:8,10:12]
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
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), '.', 'Color', color2(angle1,:),'MarkerSize',1); 
    hold on;
    scatter(axon_current(:,3),axon_current(:,1),1,[0.5,0.5,0.5],'filled');
    hold on;
    plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
        [soma_center1(2)-U(2,1)*scale1*polarity3(icell),soma_center1(2)+U(2,1)*scale1*polarity3(icell)],'color','k','LineWidth',1)
    axis off; axis image;
    hold on;
    xlim([0,600]);
    ylim([0,1200]);
end
hold on;
plot([200,400],[1100,1100],'k');
text(300,1150,'2 mm');

ax3 = subplot(1,5,4);
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,T1.bias_angle);
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = T1.soma_center;
cell_n = size(T1,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
axis off; axis image;
xlim([0,600]);
ylim([0,1200]);
scale1 = 15;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-T1.axon_bias(i,1)*scale1*T1.pc_ratio(i),soma_center1(i,1)+T1.axon_bias(i,1)*scale1*T1.pc_ratio(i)],...
        [soma_center1(i,2)-T1.axon_bias(i,2)*scale1*T1.pc_ratio(i),soma_center1(i,2)+T1.axon_bias(i,2)*scale1*T1.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
% for i = 1:cell_n
%     hold on;
%     plot([soma_center1(i,1)-T1.axon_bias(i,1)*scale1*T1.pc_ratio(i),soma_center1(i,1)+T1.axon_bias(i,1)*scale1*T1.pc_ratio(i)],...
%         [soma_center1(i,2)-T1.axon_bias(i,2)*scale1*T1.pc_ratio(i),soma_center1(i,2)+T1.axon_bias(i,2)*scale1*T1.pc_ratio(i)],...
%         'color','k','LineWidth',1);
% end
hold on;
scatter(244,542,36,'*','k');
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-pi/2','0','pi/2'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.04 a(2)+0.05 a(3) a(4)/2])% To change size

ax4 = subplot(1,5,5);
clear ang_diff1 ang_diff2
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
center = [244,542]; % SSp-un
orthog_vector = soma_center1-center;
for k = 1:100
    rp_indx = randperm(size(T1.axon_bias,1));
    vector_all1 = T1.axon_bias(rp_indx,:);
    for i =1:size(T1.axon_bias,1)
        u = [vector_all1(i,:) 0];
        v = [orthog_vector(i,:) 0];
        ang_diff1(i) = atan2(norm(cross(u,v)),dot(u,v));
    end    
    ang_diff2 = rad2deg(ang_diff1);
    edges = 0:10:180;
    [N_rp(k,:),edges] = histcounts(ang_diff2,edges);
end
ax4 = plot_angle_stairs_shuffle(ax4,edges,N1,N_rp);
hold on;
xline(mean_bias,'r');
xlim([0,180]); 
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};
% ylim([0,40]);