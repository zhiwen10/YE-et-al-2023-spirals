function h2df = plotPlaneWaveIndex(T,data_folder)
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
%%
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
    vxyRaw = complex(vxRaw, vyRaw);
    vxyRaw_abs = abs(vxyRaw);
    vxyRaw_abs = reshape(vxyRaw_abs,size(vxyRaw_abs,1),size(vxyRaw_abs,2)*size(vxyRaw_abs,3));
    vxyRaw_abs = vxyRaw_abs';
    vxyRaw_abs(not(BW1),:) = nan;
    vxyRaw_abs(vxyRaw_abs<=0.1) = nan;
    vxyRaw_abs = vxyRaw_abs';
    vxyRaw_abs = reshape(vxyRaw_abs,size(vxyRaw,1),size(vxyRaw,2),size(vxyRaw,3));
    %%
    waveIndex(kk,1) = abs(sum(vxyRaw(not(isnan(vxyRaw_abs)))))./sum(vxyRaw_abs(:),'omitnan');
end
%%
waveIndex_mean = mean(waveIndex);
waveIndex_sem = std(waveIndex)./sqrt(15);
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
