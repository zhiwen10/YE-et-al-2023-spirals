function hs10e = plotAxonFlowMatch_large2_s(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%
colora2= colorcet('C06','N',180);
%
BW = logical(projectedAtlas1);
root1 = '/997/';
ctx = '/997/8/567/688/';
scale = 8;
BW1 = BW(1:scale:end,1:scale:end);
% mask and Kernel regression map for right and left
BW1_left = BW1;
BW1_left(:,72:end) = 0;
%
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
%% get axon bias vectors for all cells with morphology
T1 = readtable(fullfile(data_folder,'axons','Axon_bias_all_cells.csv'));
soma_center_x = T1.soma_center_1;
soma_center_y = T1.soma_center_2;
axon_vector_x = T1.axon_bias_1;
axon_vector_y = T1.axon_bias_2;
pc_ratio = T1.pc_ratio;
scale = 8;
soma_center = [soma_center_x, soma_center_y]/scale;
axon_vector = [axon_vector_x, axon_vector_y];
soma_center1  = round(soma_center);
cell_n = size(soma_center,1);
for i = 1:cell_n
    axon_vector1(i,:) = axon_vector(i,:)*pc_ratio(i);
end
%% load spiral phase maps for all sessions and cancatenate
spiral_phase_all = [];
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'axons','spirals_100pixels_mean_flow',...
        [fname '_mean_flow_large.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
end
% average the phase maps and get mean spiral optical flow field
spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);   
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
vxRaw = squeeze(vxRaw); vyRaw = squeeze(vyRaw);
% sample function flow field at the locations of the cells with morphology 
vxRaw1 = vxRaw; vyRaw1 = vyRaw;
vxRaw1(:,144/2:-1:1) = -vxRaw(:,142/2+1:end);
vyRaw1(:,144/2:-1:1) = vyRaw(:,142/2+1:end);
flow_vxy = [];
for i = 1:cell_n
    flow_vxy(i,1) = vxRaw1(soma_center1(i,2),soma_center1(i,1));
    flow_vxy(i,2) = vyRaw1(soma_center1(i,2),soma_center1(i,1));
end
%%  calculate dot product 
% dot products between flow field and axon bias vector,
% for raw axon vector data and 1000x permutation of axon vectors
sum_dot = 0;
for i = 1:size(flow_vxy ,1)
    dot_a = abs(dot(flow_vxy(i,:),axon_vector1(i,:)));
    sum_dot = sum_dot+dot_a;
end
%
N_flow = vecnorm(flow_vxy,2,2);
N_axon = vecnorm(axon_vector1,2,2);
N_abs = sum(N_flow.*N_axon);
%
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
[h1,p] = ttest2(sum_dot,sum_dot_perm_all);
%% plot optical flow fields and matching index
flow_angle = round(rad2deg(atan(vyRaw1(:)./vxRaw1(:))));
x1 = linspace(-90,90,180);
colora3b = interp1(x1,colora2,flow_angle);
colora3b = reshape(colora3b,[165,143,3]);
vxRaw2 = nan(size(vxRaw1));
vyRaw2 = nan(size(vyRaw1));
skip = 3; zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw1(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw1(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
[row,col] = find(not(isnan(vyRaw2)));
row = row(col<71);col = col(col<71); % only plot left column
lineColor = 'k'; lineColor1 = 'w';
hemi = 'left';
scale3 = 5/8;
%%
center = [350,700]; % PPC
hs10e = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
ax1 = subplot(1,3,1);
x1 = linspace(-90,90,180);
color2= colorcet('C06','N',180);
colora3 = interp1(x1,color2,T1.bias_angle);
colora3(isnan(colora3)) = 0;
vector_scale = 2;ring_scale = 10;scale = 1; 
soma_center2 = [T1.soma_center_1,T1.soma_center_2];
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
scale1 = 2;
for i = 1:cell_n
% for i = 1:cell_n
    hold on;
    plot([soma_center1(i,1)-axon_bias(i,1)*scale1*T1.pc_ratio(i),...
        soma_center1(i,1)+axon_bias(i,1)*scale1*T1.pc_ratio(i)],...
        [soma_center1(i,2)-axon_bias(i,2)*scale1*T1.pc_ratio(i),...
        soma_center1(i,2)+axon_bias(i,2)*scale1*T1.pc_ratio(i)],...
        'color',colora3(i,:),'LineWidth',1);
end
hold on;
scatter(center(1)/8,center(2)/8,36,'*','k');

ax2 = subplot(1,3,2);
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
soma_n = size(soma_center1,1);
scale1 = 1;
hold on;
for i = 1:soma_n
    plot([soma_center1(i,1)-flow_vxy(i,1)*scale1,soma_center1(i,1)+flow_vxy(i,1)*scale1],...
        [soma_center1(i,2)-flow_vxy(i,2)*scale1,soma_center1(i,2)+flow_vxy(i,2)*scale1],...
        'color',squeeze(colora3b(soma_center1(i,2),soma_center1(i,1),:)),'LineWidth',1);
    hold on;
end
hold on;
scatter(center(1)/8,center(2)/8,36,'*','k');
set(gca,'Ydir','reverse')
axis image; axis off;

sum_dot_perm_all1 = sum_dot_perm_all./N_abs;
sum_dot1 = sum_dot/N_abs;
[h2,p2] = ttest2(sum_dot1,sum_dot_perm_all1);
%

% [angle_diff_real,angle_diff_perm,edges,N1,N_rp] = getAngleHistogramStats(T1);
% [h mu ul ll] = circ_mtest(angle_diff_real,deg2rad(90));
% mean_rad = [mu,ll,ul];
% mean_deg = rad2deg(mean_rad);
% for i = 1:100
%     [p(i)]=watsons_U2_perm_test(angle_diff_real,angle_diff_perm(i,:),100); % good
% end
% p = mean(p);

ax3 = subplot(1,3,3);
histogram(sum_dot_perm_all1,'faceColor','None','EdgeColor','k');
hold on;
xline(sum_dot1);
text(1,120,num2str(p2));
text(1,130,num2str(sum_dot1));
ylim([0,200]);
xlim([0.5,0.8]);
% xticks([0.5:0.1:0.8]);
% xticklabels({'0.5','0.6','0.7','0.8'});
yticks([0:50:200]);
yticklabels({'0','50','100','150','200'});

print(hs10e, fullfile(save_folder,'FigS10e_axon_flow_match_large2.pdf'),...
    '-dpdf', '-bestfit', '-painters');