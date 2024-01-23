T1 = readtable('axon_bias_table2.csv');
% a1 = not(ismember(T1.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"])));
% T1 = T1(a1,:);
soma_center_x = T1.soma_center_1;
soma_center_y = T1.soma_center_2;
axon_vector_x = T1.axon_bias_1;
axon_vector_y = T1.axon_bias_2;
pc_ratio = T1.pc_ratio;
%%
% load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralExamples\sprial_frame_optical_flow.mat');
% load('C:\Users\Steinmetz lab\Documents\git\BIL\function_matching\average_flow_phase_mean.mat');
load('C:\Users\Steinmetz lab\Documents\git\BIL\function_matching\average_flow2.mat');
BW = BW1; vxRaw = vxRaw_mean; vyRaw = vyRaw_mean; 
%%
vxRaw1 = vxRaw; vyRaw1 = vyRaw;
vxRaw1(:,144/2:-1:1) = -vxRaw(:,142/2+1:end);
vyRaw1(:,144/2:-1:1) = vyRaw(:,142/2+1:end);
%%
flow_angle = round(rad2deg(atan(vyRaw_mean(:)./vxRaw_mean(:))));
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,flow_angle);
colora3 = reshape(colora3,[165,143,3]);
%%
scale = 8;
soma_center = [soma_center_x, soma_center_y]/scale;
axon_vector = [axon_vector_x, axon_vector_y];
%%
figure;
ax1 = subplot(1,3,1);
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
skip = 3; zoom_scale = 5;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW) = nan; vyRaw2(~BW) = nan;
scale2 = 8;
overlayOutlines(coords,scale2);
scale1 = 0.5;
for i = 1:numel(row)
    plot([col(i)-vxRaw2(row(i),col(i))*scale1,col(i)+vxRaw2(row(i),col(i))*scale1],...
        [row(i)-vyRaw2(row(i),col(i))*scale1,row(i)+vyRaw2(row(i),col(i))*scale1],...
        'color',colora3(row(i),col(i),:),'LineWidth',1);
    hold on;
end
axis image; axis off;
set(gca,'Ydir','reverse')
%
ax2 = subplot(1,3,2);
scale2 = 8;
overlayOutlines(coords,scale2);
soma_center1 = round(soma_center);
scale1 = 3;
cell_n = size(soma_center,1);
for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_vector(i,1)*scale1*pc_ratio(i),soma_center1(i,1)+axon_vector(i,1)*scale1*pc_ratio(i)],...
        [soma_center1(i,2)-axon_vector(i,2)*scale1*pc_ratio(i),soma_center1(i,2)+axon_vector(i,2)*scale1*pc_ratio(i)],...
        'color','k','LineWidth',0.5);
end
axis image; axis off;
set(gca,'Ydir','reverse')
%
ax3 = subplot(1,3,3);
scale2 = 8;
overlayOutlines(coords,scale2);
axis image; axis off;
soma_mirror = soma_center1;
soma_mirror(:,1) = 143-soma_mirror(:,1);
flow_vxy = [];
for i = 1:cell_n
    flow_vxy_plot(i,1) = vxRaw1(soma_mirror(i,2),soma_mirror(i,1));
    flow_vxy_plot(i,2) = vyRaw1(soma_mirror(i,2),soma_mirror(i,1));
end
scale1 = 8;
cell_n = size(soma_center,1);
for i = 1:cell_n
    hold on;
    quiver(soma_mirror(:,1),soma_mirror(:,2),flow_vxy_plot(:,1)*scale1,flow_vxy_plot(:,2)*scale1,'k','lineWidth',0.5,'autoScale','off');
end
set(gca,'Ydir','reverse')
%%
for i = 1:cell_n
    axon_vector1(i,:) = axon_vector(i,:)*pc_ratio(i);
end
flow_vxy = [];
for i = 1:cell_n
    flow_vxy(i,1) = vxRaw1(soma_center1(i,2),soma_center1(i,1));
    flow_vxy(i,2) = vyRaw1(soma_center1(i,2),soma_center1(i,1));
end
sum_dot = 0;
for i = 1:size(flow_vxy ,1)
    dot_a = abs(dot(flow_vxy(i,:),axon_vector1(i,:)));
    sum_dot = sum_dot+dot_a;
end
sum_dot_perm_all = [];
for k = 1:1000
    index = randperm(size(flow_vxy ,1));
    axon_vector_perm = axon_vector1(index,:);
    sum_dot_perm = 0;
    for i = 1:size(flow_vxy ,1)
        dot_a = abs(dot(flow_vxy(i,:),axon_vector_perm(i,:)));
        sum_dot_perm = sum_dot_perm+dot_a;
    end
    sum_dot_perm_all(k,1) = sum_dot_perm;
end
figure;
histogram(sum_dot_perm_all);
hold on;
xline(sum_dot);
[h1,p] = ttest2(sum_dot,sum_dot_perm_all);
%% mask and Kernel regression map for right and left
BW1_left = BW1;
BW1_left(:,72:end) = 0;
%% cortex surface outline
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
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
%%
h1 = figure;
framea1 = fliplr(framea);
ax1 = subplot(1,2,1);
hold off;
im_phase = imagesc(framea1);
colormap(ax1,colorcet('C06'));
axis image; axis off;
scale1 = 8;
hold on;
quiver(soma_center1(:,1),soma_center1(:,2),flow_vxy(:,1)*scale1,flow_vxy(:,2)*scale1,'k','lineWidth',0.5,'autoScale','off');
set(im_phase, 'AlphaData', BW1_left, 'AlphaDataMapping', 'scaled');
lineColor = 'w'; lineColor1 = 'w';
hemi = 'left';
scale3 = 5/8;
hold on;
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);

ax2 = subplot(1,2,2)
histogram(sum_dot_perm_all);
hold on;
xline(sum_dot);
[h1,p] = ttest2(sum_dot,sum_dot_perm_all);
text(140,120,num2str(p));
% print(h1, 'mean_flow_vector', '-dpdf', '-bestfit', '-painters');