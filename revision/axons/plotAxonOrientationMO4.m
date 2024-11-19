function [hr3ad] = plotAxonOrientationMO4(data_folder,save_folder)
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
load(fullfile(data_folder,'Revision','axons','MO_roi.mat'));
%%
center = [377,428]; % MOp
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
hr3ad = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
ax1 = subplot(1,4,1);
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO.csv'));
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
scatter(center(1),center(2),36,'*','k');
hold on;
scatter(center2(1),center2(2),36,'*','k');

ax2 = subplot(1,4,2);
T1 = readtable(fullfile(data_folder,'Revision','axons','Axon_bias_all_cells_MO2.csv'));
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
scatter(center(1),center(2),36,'*','k');
hold on;
plot([point1(1),point2(1)],[point1(2),point2(2)],'--k','LineWidth',1);
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
ylim([0,15]);
xticks =  [0,45,90,135,180];
xticklabels = {'0','1/4*pi','1/2*pi','3/4*pi','pi'};
% define border by point1 and point2
slope = (point1(2)-point2(2))./(point1(1)-point2(1));
intercept = point1(2)-slope*point1(1);
x = T1.soma_center_1;
y = T1.soma_center_2;
border = y- (x*slope+intercept);
left_roi = (border<0);
right_roi = (border>=0);
% T4 = T1(left_roi,:);
center = [377,428]; % MOp
axon_vector = [T1.axon_bias_1,T1.axon_bias_2];
theta = 90;
for i = 1:size(axon_vector,1)
    v1 = axon_vector(i,:);
    vr = v1*[1;1i]*exp(-1i*theta*pi/180);
    angle1(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
end
angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
angle2a = angle1(left_roi);
angle2b = angle1(right_roi);
%
ax4 =subplot(1,4,4);
hist2a = histogram(angle2a,10,'faceColor','None','EdgeColor','k','DisplayStyle','stairs');
angle2a_max = mean(hist2a.BinEdges(hist2a.Values == max(hist2a.Values)));
xline(angle2a_max+5,'k--');
hold on;
hist2b = histogram(angle2b,10,'faceColor','None','EdgeColor','r','DisplayStyle','stairs');
angle2b_max = hist2b.BinEdges(hist2b.Values == max(hist2b.Values));
xline(angle2b_max+5,'k--');
xlim([-90,90]);
ylim([0,30]);
ax4.XTick =  [-90,-45,0,45,90];
ax4.XTickLabel = {'-pi/2','-pi/4','0','pi/4','pi/2'};
%
angle_all = [angle2a_max+6,angle2b_max+6];
angle_all2 = angle_all+90;
angle_all2 = angle_all2-180;
color2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,color2,angle_all2);
hist2a.EdgeColor = colora3(1,:);
hist2b.EdgeColor = colora3(2,:);
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
print(hr3ad, fullfile(save_folder,'FigR3ad_circularbias_MO2.pdf'),...
    '-dpdf', '-bestfit', '-painters');