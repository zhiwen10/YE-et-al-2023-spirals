%% get atlas mask and outlines
% load coords for atlas
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
% load projectedAtlas and projectedTemplate
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
BW = logical(projectedAtlas1);
scale = 8;
BW1 = BW(1:scale:end,1:scale:end);
%% mask and Kernel regression map for right and left
BW1_left = BW1;
BW1_left(:,72:end) = 0;
%% cortex surface outline
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
root1 = '/997/';
areaName = {'MOp','MOs'};
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
T1 = readtable('C:\Users\Steinmetz lab\Documents\git\BIL\axon_bias_table2.csv');
% a1 = not(ismember(T1.labels,string(["VIS";"RSP";"AUD";"TEa";"VISC"])));
% T1 = T1(a1,:);
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
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
%%
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
data_folder = 'C:\Users\Steinmetz lab\Documents\git\BIL\function_matching\spirals_70pixels';
spiral_phase_all = [];
vxRaw_all1 = [];
vyRaw_all1 = [];
for kk = 1:15
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,[fname '_mean_phase_flow.mat']),'spiral_phase_all_norm','vxRaw_all','vyRaw_all');
    spiral_phase_all = cat(3,spiral_phase_all,spiral_phase_all_norm);
    vxRaw_all1 = cat(3,vxRaw_all1,vxRaw_all);
    vyRaw_all1 = cat(3,vyRaw_all1,vyRaw_all);
end
%%
spiral_phase_mean = circ_mean(spiral_phase_all, [], 3);
vxRaw = squeeze(mean(vxRaw_all,3)); vyRaw = squeeze(mean(vyRaw_all,3));
%%
vxRaw1 = vxRaw; vyRaw1 = vyRaw;
vxRaw1(:,144/2:-1:1) = -vxRaw(:,142/2+1:end);
vyRaw1(:,144/2:-1:1) = vyRaw(:,142/2+1:end);
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
[h1,p] = ttest2(sum_dot,sum_dot_perm_all);
%%
flow_angle = round(rad2deg(atan(vyRaw1(:)./vxRaw1(:))));
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,flow_angle);
colora3 = reshape(colora3,[165,143,3]);
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
h3 = figure;
ax1 = subplot(1,3,1);
for i = [1:15,19:22]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor1);
end
plotOutline(maskPath([1:6,14,15,19]),st,atlas1,hemi,scale3,lineColor);
plotOutline(maskPath([16:18]),st,atlas1,hemi,scale3,lineColor);
for i = [14:15]
    plotOutline(maskPath(i),st,atlas1,hemi,scale3,lineColor);
end
hold on;
scale1 = 1;
for i = 1:numel(row)
    plot([col(i)-vxRaw2(row(i),col(i))*scale1,col(i)+vxRaw2(row(i),col(i))*scale1],...
        [row(i)-vyRaw2(row(i),col(i))*scale1,row(i)+vyRaw2(row(i),col(i))*scale1],...
        'color',colora3(row(i),col(i),:),'LineWidth',1);
    hold on;
end
hold on;
scatter(244/8,542/8,36,'*','k');
set(gca,'Ydir','reverse')
axis image; axis off;

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
scale1 = 2;
hold on;
for i = 1:soma_n
    plot([soma_center1(i,1)-flow_vxy(i,1)*scale1,soma_center1(i,1)+flow_vxy(i,1)*scale1],...
        [soma_center1(i,2)-flow_vxy(i,2)*scale1,soma_center1(i,2)+flow_vxy(i,2)*scale1],...
        'color',colora3(soma_center1(i,2),soma_center1(i,1),:),'LineWidth',1);
    hold on;
end
% quiver(soma_center1(:,1),soma_center1(:,2),flow_vxy(:,1)*scale1,flow_vxy(:,2)*scale1,'k','lineWidth',0.5,'autoScale','off');
hold on;
scatter(244/8,542/8,36,'*','k');
set(gca,'Ydir','reverse')
axis image; axis off;

min_indx = min(sum_dot_perm_all);
sum_dot_perm_all1  = (sum_dot_perm_all-min_indx);
max_indx = max(sum_dot_perm_all1);
sum_dot_perm_all1  = sum_dot_perm_all1/max_indx;
sum_dot1 = (sum_dot-min_indx)/max_indx;
[h2,p2] = ttest2(sum_dot1,sum_dot_perm_all1);
ax2 = subplot(1,3,3)
histogram(sum_dot_perm_all1);
hold on;
xline(sum_dot1);
text(1,120,num2str(p2));
text(1,130,num2str(sum_dot1));
ylim([0,150]);
print(h3, ['matching_index_across_sessions'], '-dpdf', '-bestfit', '-painters');