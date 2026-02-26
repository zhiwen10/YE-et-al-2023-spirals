function [h2ac] = plotAxonOrientation2(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells2.csv'));
load(fullfile(data_folder, 'axons','example_12_cells.mat'));
%% examples
h2ac = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
%
subplot(1,4,1);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
hold on;
center = [244,542]; % SSp-un
% center = [236,629]; % SSp-un
scatter(center(1), center(2),36,'*','k');
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

subplot(1,4,2);
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for icell = 1
    axon_current_all = axon_current_all1{icell};
    axon_current = axon_current1{icell};
    soma_current_2d = soma_current_2d1{icell};
    
    axon_current_2d = double(axon_current(:,[3,1]));
    
    axon_center = mean(axon_current_2d,1);
    axon_current_2d = axon_current_2d-...
        repmat(axon_center,size(axon_current_2d,1),1);    
    soma_center1 = double(soma_current_2d(:,[3,1]));
    sample_n = size(axon_current_2d,1);
    [U,S,V] = svd(double(axon_current_2d)',"econ");
    S_diag3 = diag(S);
    angle1 = round(rad2deg(atan(U(2,1)/U(1,1))))+91;
    polarity3(icell) = S_diag3(1,:)./S_diag3(2,:)-1;
    hold on;
    plot(axon_current_all(:,3), axon_current_all(:,1), ...
        '.', 'Color', color2(angle1,:),'MarkerSize',1); 
    hold on;
    scatter(axon_current(:,3),axon_current(:,1),1,[0.5,0.5,0.5],'filled');
    hold on;
    plot([soma_center1(1)-U(1,1)*S_diag3(1)/sqrt(sample_n),...
        soma_center1(1)+U(1,1)*S_diag3(1)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,1)*S_diag3(1)/sqrt(sample_n),...
        soma_center1(2)+U(2,1)*S_diag3(1)/sqrt(sample_n)],...
        'color','k','LineWidth',1);
    hold on;
    plot([soma_center1(1)-U(1,2)*S_diag3(2)/sqrt(sample_n),...
        soma_center1(1)+U(1,2)*S_diag3(2)/sqrt(sample_n)],...
        [soma_center1(2)-U(2,2)*S_diag3(2)/sqrt(sample_n),...
        soma_center1(2)+U(2,2)*S_diag3(2)/sqrt(sample_n)],...
        'color','k','LineWidth',1);
    axis off; 
    axis image;
end
hold on;
scatter(center(1), center(2),36,'*','k');
hold on;
plot([soma_center1(1),center(1)],[soma_center1(2),center(2)],'--','color','k');
hold on;
plot([150,250],[800,800],'k');
text(150,850,'1 mm');
xlim([100,500]);
ylim([400,800]);

subplot(1,4,3);
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

ax3 = subplot(1,4,4);
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
% scale1 = 15;
scale1 = 0.1;
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*T1.axon_length(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*T1.axon_length(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*T1.axon_length(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*T1.axon_length(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1), center(2),36,'*','k');
cb = colorbar;
colormap(flipud(color2));
cb.Ticks =  [0,0.5,1];
cb.TickLabels = {'-pi/2','0','pi/2'};
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1)+0.04 a(2)+0.05 a(3) a(4)/2])% To change size
%%
print(h2ac, fullfile(save_folder,'Fig2ac_axon_bias_all.pdf'),...
    '-dpdf', '-bestfit', '-painters');