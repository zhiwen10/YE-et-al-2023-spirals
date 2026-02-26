function [hs10ad] = plotAxonOrientationLayer(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
% mask and Kernel regression map for right and left
BW = logical(projectedAtlas1);
BW1_left = BW;
BW1_left(:,72:end) = 0;
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
%%
center2 = [244,542]; % SSp-un
point1 = [340,475];
point2 = [500,250];
scale = 1;
scale1 = 15;
color1= colorcet('C06','N',12);
color2= colorcet('C06','N',180);
%%
scale3 = 5;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'left';
hs10ad = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
ax1 = subplot(1,4,1);
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells3.csv'));
T1 = T1(any(T1{:,11:13},2),:);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
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
scatter(center2(1),center2(2),36,'*','k');

ax4 =subplot(1,4,3);
clear ang_diff1 ang_diff2
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
center = [244,542]; % SSp-un
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
ylim([0,25]);
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};

[angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStats(T1);
[h1 mu1 ul1 ll1] = circ_mtest(angle_diff_real,deg2rad(90));
mean_rad1 = [mu1,ll1,ul1];
mean_deg1 = rad2deg(mean_rad1);
clear p1
for i = 1
    [p1(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
end

ax2 = subplot(1,4,2);
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells3.csv'));
T1 = T1(any(T1{:,14:15},2),:);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,T1.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
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
scatter(center2(1),center2(2),36,'*','k');

ax4 =subplot(1,4,4);
clear ang_diff1 ang_diff2
soma_center1 = [T1.soma_center_1,T1.soma_center_2];
axon_bias = [T1.axon_bias_1,T1.axon_bias_2];
edges = 0:10:180;
[N1,edges] = histcounts(T1.center_bias_angle,edges);   
mean_bias = mean(T1.center_bias_angle);
% shuffle
center = [244,542]; % SSp-un
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
ylim([0,25]);
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};

[angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStats(T1);
[h2 mu2 ul2 ll2] = circ_mtest(angle_diff_real,deg2rad(90));
mean_rad2 = [mu2,ll2,ul2];
mean_deg2 = rad2deg(mean_rad2);
clear p2
i = 1;
[p2(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
%%
print(hs10ad, fullfile(save_folder,'FigS10ad_layers.pdf'),...
    '-dpdf', '-bestfit', '-painters');