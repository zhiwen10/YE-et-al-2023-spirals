function h2df = plotFlowFieldSummary(T,data_folder,save_folder)
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
%% load spiral phase maps for all sessions and cancatenate
spiral_phase_all = [];
for kk = 1:15
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'axons','spirals_70pixels_mean_flow',...
        [fname '_mean_flow2.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
    load(fullfile(data_folder,'axons','spirals_70pixels_mean_flow_left',...
        [fname '_mean_flow2_left.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
end
% average the phase maps and get mean spiral optical flow field
spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);   
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
vxRaw = squeeze(vxRaw); vyRaw = squeeze(vyRaw);
%% load spiral phase maps for all sessions and cancatenate
vxyRaw_abs = nan(165,143,15);
for kk = 1:size(T,1)
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,[fname '_spiral_flow.mat']),'vxRaw_mean','vyRaw_mean');
    %%
    vxyRaw = complex(vxRaw_mean,vyRaw_mean);
    vxyRaw_abs(:,:,kk) = abs(vxyRaw);
end
%%
pixel(1,:) = [97,86]; % RSP
pixel(2,:) = [110,110]; % VISp
pixel(3,:) = [60,80]; % ACC
pixel(4,:) = [85,115]; % SSp

h1 = figure('Renderer', 'painters', 'Position', [100 100 1200 400]);
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
hold on;
scatter(pixel(:,2),pixel(:,1));

ax2 = subplottight(1,4,2);
frameb = squeeze(spiral_phase_mean(:,:,2));
im_phase = imagesc(frameb);
colormap(ax2,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
scatter(pixel(:,2),pixel(:,1));

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
hold on;
scatter(pixel(:,2),pixel(:,1));

ax4 = subplottight(1,4,4);
vxyRaw = complex(vxRaw,vyRaw);
vxyRaw_abs = abs(vxyRaw);
im_phase = imagesc(vxyRaw_abs);
colormap(ax4,parula);
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
scatter(pixel(:,2),pixel(:,1));

print(h1,'allspirals_average.pdf','-dpdf', '-bestfit', '-painters');