function h2df = plotFlowFieldSummary2(T,data_folder,save_folder)
%%
save_folder = 'E:\spiral_data_share\data\revision2';
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
BW = logical(projectedAtlas1);
root1 = '/997/';
ctx = '/997/8/567/688/';
scale = 8;
BW1 = BW(1:scale:end,1:scale:end);
%% get speed and flow vector length scale
kk = 15;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
% load raw detected spirals
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(save_folder,'right',[fname '_spiral_flow.mat']));
vxyRaw = complex(vxRaw_mean,vyRaw_mean);
vxyRaw_abs1 = abs(vxyRaw);
phase_mean_all = squeeze(phase_mean_all);
pixel_ssp = [85,115];
phase_example = squeeze(phase_mean_all(:,pixel_ssp(1),pixel_ssp(2)));
phase_diff = phase_example(2)-phase_example(1);
spiral_center = [72,110];
angular_speed = phase_diff*35;
point_diff = pixel_ssp-spiral_center;
radius = vecnorm(point_diff,2,2)*8*0.01;
linear_speed = angular_speed*radius;
vector_length = squeeze(vxyRaw_abs1(pixel_ssp(1),pixel_ssp(2)));
speed_scale = linear_speed./vector_length; %mm/s
%%
figure;
ax1 = subplot(1,2,1);
framea = squeeze(phase_mean_all(1,:,:));
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; 
ax2 = subplot(1,2,2);
framea = squeeze(phase_mean_all(2,:,:));
im_phase = imagesc(framea);
colormap(ax2,colorcet('C06'));
axis image; 
%%
pixel(1,:) = [97,86]; % RSP
pixel(2,:) = [110,110]; % VISp
pixel(3,:) = [60,80]; % ACC
pixel(4,:) = [85,115]; % SSp
%%
vector_norm = [];
for kk = 1:size(T,1)
    %%
    clear vector_all
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,'right',[fname '_spiral_flow.mat']));
    %%
    for i = 1:4
        vector_all(i,1,:) = squeeze(vxRaw(:,pixel(i,1),pixel(i,2)));
        vector_all(i,2,:) = squeeze(vyRaw(:,pixel(i,1),pixel(i,2)));
        vector2(kk,i,1) = vxRaw_mean(pixel(i,1),pixel(i,2));
        vector2(kk,i,2) = vyRaw_mean(pixel(i,1),pixel(i,2));
    end
    %%
    vector_all_length = vecnorm(vector_all,2,2);
    vector_all_norm = vector_all./vector_all_length;
    vector_norm(kk,:) = vecnorm(sum(vector_all_norm,3),2,2)./size(vector_all_norm,3);
end
%%
vector2(:,:,2) = -vector2(:,:,2);
vector2 = vector2.*speed_scale;
speed = vecnorm(vector2,2,3);
%%
color1 = {'r','b','g','m'};
h1= figure('Renderer', 'painters', 'Position', [100 100 800 200]);
subplot(1,3,1);
x = zeros(15,1); y = zeros(15,1);
compass(squeeze(vector2(:,4,1)),squeeze(vector2(:,4,2)),color1{4});
hold on;
compass(squeeze(vector2(:,1,1)),squeeze(vector2(:,1,2)),color1{1});
hold on;
compass(squeeze(vector2(:,2,1)),squeeze(vector2(:,2,2)),color1{2});
hold on;
compass(squeeze(vector2(:,3,1)),squeeze(vector2(:,3,2)),color1{3});

subplot(1,3,2);
pairs = [1,2;1,3;1,4;2,3;2,4;3,4];
for pixel_i = 1:4
    vector_i = squeeze(speed(:,pixel_i));    
    id_i = ones(size(vector_i))*pixel_i;     
    scatter(id_i,vector_i,4,color1{pixel_i},'filled');
    hold on;
end
pairs2 = [1,2;2,3;3,4];
for i = 1:size(pairs2,1)
    vector_i1 = squeeze(speed(:,pairs2(i,1)));  
    vector_i2 = squeeze(speed(:,pairs2(i,2)));  
    for kk = 1:15         
        plot(pairs2(i,:)',[vector_i1(kk),vector_i2(kk)]','k');
        hold on;
    end
end
xticks([1,2,3,4]);
xticklabels({'RSP','VISp','ACC','SSp'});
yticks([0:10:50]);
yticklabels([0:10:50]);
ylim([0,50]);
ylabel('Linear speed (mm/s)');

subplot(1,3,3);
pairs = [1,2;1,3;1,4;2,3;2,4;3,4];
for pixel_i = 1:4
    vector_i = squeeze(vector_norm(:,pixel_i));    
    id_i = ones(size(vector_i))*pixel_i;     
    scatter(id_i,vector_i,4,color1{pixel_i},'filled');
    hold on;
end
pairs2 = [1,2;2,3;3,4];
for i = 1:size(pairs2,1)
    vector_i1 = squeeze(vector_norm(:,pairs2(i,1)));  
    vector_i2 = squeeze(vector_norm(:,pairs2(i,2)));  
    for kk = 1:15         
        plot(pairs2(i,:)',[vector_i1(kk),vector_i2(kk)]','k');
        hold on;
    end
end
xticks([1,2,3,4]);
xticklabels({'RSP','VISp','ACC','SSp'});
yticks([0:0.2:1]);
yticklabels([0:0.2:1]);
ylim([0,1]);
ylabel('Coherence');
print(h1, fullfile(save_folder,'flow_vector_speed_cohenrence.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%%
figure;
hemi = [];
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';

vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
skip = 5;
zoom_scale = 2;
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;

ax1 = subplottight(1,4,1);
framea = squeeze(spiral_phase_mean(:,:,1));
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');

ax2 = subplottight(1,4,2);
frameb = squeeze(spiral_phase_mean(:,:,2));
im_phase = imagesc(frameb);
colormap(ax2,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');

ax3 = subplottight(1,4,3);
framea = squeeze(spiral_phase_mean(:,:,1));
im_phase = imagesc(framea);
colormap(ax3,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');

ax4 = subplottight(1,4,4);
vxyRaw = complex(vxRaw,vyRaw);
vxyRaw_abs = abs(vxyRaw);
im_phase = imagesc(vxyRaw_abs);
colormap(ax4,parula);
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
