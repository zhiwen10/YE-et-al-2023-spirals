function h2df = plotFlowFieldAcrossSessions(T,data_folder,save_folder)
%%
save_folder = 'E:\spiral_data_share\data\revision2\left';
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
useGPU = 0;
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'axons','spirals_70pixels_mean_flow_left',...
        [fname '_mean_flow2_left.mat']),'spiral_phase_all_norm');
    spiral_phase_all_norm = permute(spiral_phase_all_norm,[4,3,1,2]);
    %%
    vxRaw = zeros(size(spiral_phase_all_norm,1),...
        size(spiral_phase_all_norm,3),size(spiral_phase_all_norm,4));
    vyRaw = zeros(size(spiral_phase_all_norm,1),...
        size(spiral_phase_all_norm,3),size(spiral_phase_all_norm,4));
    for i = 1:size(spiral_phase_all_norm,1)
        phase_temp = squeeze(spiral_phase_all_norm(i,:,:,:));
        A1 = squeeze(phase_temp(1,:,:));
        A2 = squeeze(phase_temp(2,:,:));
        [vxRaw(i,:,:), vyRaw(i,:,:)] = HS_phase_mod(A1, A2);
    end
    vxRaw_mean = squeeze(mean(vxRaw,1));
    vyRaw_mean = squeeze(mean(vyRaw,1));
    phase_mean = circ_mean(squeeze(spiral_phase_all_norm(:,1,:,:)),[],1);
    phase_mean_all = circ_mean(spiral_phase_all_norm,[],1);
    %%
    save(fullfile(save_folder,[fname '_spiral_flow_left.mat']),'vxRaw','vyRaw',...
        'vxRaw_mean','vyRaw_mean','phase_mean_all');
end
%%
figure;
hemi = [];
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';

ax1 = subplottight(1,1,1);
framea = squeeze(phase_mean);
vxRaw2 = nan(size(vxRaw_mean));
vyRaw2 = nan(size(vyRaw_mean));
skip = 5;
zoom_scale = 10;
vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;

im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');