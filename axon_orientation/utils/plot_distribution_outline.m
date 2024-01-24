function h = plot_distribution(T1,coords,st)
%% cortex surface outline
% st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
root1 = '/997/';
maskPath{1} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{2} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{3} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{4} = '/997/8/567/688/695/315/22'; % VISa
maskPath{5} = '/997/8/567/688/695/315/677/'; % VISC
maskPath{6} = '/997/8/567/688/695/315/22/417/'; % VISrl
maskPath{7} = '/997/8/567/688/695/315/453/322/361/'; % SSp-tr
maskPath{8} = '/997/8/567/688/695/315/453/322/337/'; % SSp-ll
maskPath{9} = '/997/8/567/688/695/315/453/322/369/'; % SSp-ul
maskPath{10} = '/997/8/567/688/695/315/453/322/345/'; % SSp-m
maskPath{11} = '/997/8/567/688/695/315/453/322/353/'; % SSp-n
maskPath{12} = '/997/8/567/688/695/315/453/322/329/'; % SSp-bfd
maskPath{13} = '/997/8/567/688/695/315/453/322/182305689/'; % SSp-un
maskPath{14} = '/997/8/567/688/695/315/453/378/'; % SSs

maskPath{15} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{16} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{17} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{18} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{19} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{20} = '/997/8/567/688/695/315/669/385/'; %VISp
maskPath{21} = '/997/8/567/688/695/315/669/312782628/';
maskPath{22} = '/997/8/567/688/695/315/669/409/';
ctx = '/997/8/567/688/';
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
%% examples
h = figure('Renderer', 'painters', 'Position', [100 100 1100 300]);
load('example_12_cells.mat');
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
lineColor = 'k'; lineColor1 = 'w';
hemi = 'left';
scale3 = 5;

subplot(1,5,1);
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
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
    scalev = 3;
%     plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
%         [soma_center1(2)-U(2,1)*scale1*polarity3(icell),soma_center1(2)+U(2,1)*scale1*polarity3(icell)],'color','k','LineWidth',1);
    plot([soma_center1(1)-U(1,1)*S_diag3(1)/sqrt(sample_n)*scalev,soma_center1(1)+U(1,1)*S_diag3(1)/sqrt(sample_n)*scalev],...
        [soma_center1(2)-U(2,1)*S_diag3(1)/sqrt(sample_n)*scalev,soma_center1(2)+U(2,1)*S_diag3(1)/sqrt(sample_n)*scalev],'color','k','LineWidth',1);
    hold on;
    plot([soma_center1(1)-U(1,2)*S_diag3(2)/sqrt(sample_n)*scalev,soma_center1(1)+U(1,2)*S_diag3(2)/sqrt(sample_n)*scalev],...
        [soma_center1(2)-U(2,2)*S_diag3(2)/sqrt(sample_n)*scalev,soma_center1(2)+U(2,2)*S_diag3(2)/sqrt(sample_n)*scalev],'color','k','LineWidth',1);
    axis off; 
    axis image;
end
set(gca, 'YDir','reverse');
hold on;
center = [244,542]; % SSp-un
scatter(244,542,36,'*','k');
hold on;
plot([soma_center1(1),244],[soma_center1(2),542],'-','color','k');
hold on;
scatter(soma_center1(1),soma_center1(2),5,'r','filled');
hold on;
plot([200,400],[1100,1100],'k'); % scale bar
text(300,1150,'2 mm');
hold on;
v = [150,400;150,750;450,750;450,400]/scale;
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
axis off; 
axis image;
xlim([0,600]);
ylim([0,1200]);
%
subplot(1,5,2);
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
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
    scatter(axon_current(:,3),axon_current(:,1),3,[0.5,0.5,0.5],'filled');
    hold on;
%     plot([soma_center1(1)-U(1,1)*scale1*polarity3(icell),soma_center1(1)+U(1,1)*scale1*polarity3(icell)],...
%         [soma_center1(2)-U(2,1)*scale1*polarity3(icell),soma_center1(2)+U(2,1)*scale1*polarity3(icell)],'color','k','LineWidth',1);
    plot([soma_center1(1)-U(1,1)*S_diag3(1)/sqrt(sample_n)*scalev,soma_center1(1)+U(1,1)*S_diag3(1)/sqrt(sample_n)*scalev],...
        [soma_center1(2)-U(2,1)*S_diag3(1)/sqrt(sample_n)*scalev,soma_center1(2)+U(2,1)*S_diag3(1)/sqrt(sample_n)*scalev],'color','k','LineWidth',1);
    hold on;
    plot([soma_center1(1)-U(1,2)*S_diag3(2)/sqrt(sample_n)*scalev,soma_center1(1)+U(1,2)*S_diag3(2)/sqrt(sample_n)*scalev],...
        [soma_center1(2)-U(2,2)*S_diag3(2)/sqrt(sample_n)*scalev,soma_center1(2)+U(2,2)*S_diag3(2)/sqrt(sample_n)*scalev],'color','k','LineWidth',1);
    axis off; 
    axis image;
end
hold on;
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
hold on;
center = [244,542]; % SSp-un
scatter(244,542,36,'*','k');
hold on;
plot([soma_center1(1),244],[soma_center1(2),542],'-','color','k');
hold on;
scatter(soma_center1(1),soma_center1(2),5,'r','filled');
hold on;
plot([150,250],[720,720],'k');
text(150,800,'1 mm');
xlim([150,450]);
ylim([400,750]);

subplot(1,5,3);
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
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
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
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
end