%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
data_folder = 'C:\Users\Steinmetz lab\Documents\git\BIL\function_matching\spirals_all_mean_flow';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
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
    %%
    load(fullfile(data_folder,[fname '_mean_flow2.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
end
%% get atlas mask and outlines
% load coords for atlas
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
% load projectedAtlas and projectedTemplate
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
BW = logical(projectedAtlas1);
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
%%
spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
%%
spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);    
[vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
%%
skip = 3; zoom_scale = 2;
vxRaw_mean = squeeze(vxRaw); vyRaw_mean = squeeze(vyRaw);
vxRaw2 = nan(size(BW1));
vyRaw2 = nan(size(BW1));
vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
scale2 = 8;
h = figure;
subplot(1,2,1);
im = imagesc(squeeze(spiral_phase_mean(:,:,1)));
colormap(colorcet('C06'));
axis image; axis off;
set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
subplot(1,2,2);
im = imagesc(squeeze(spiral_phase_mean(:,:,2)));
colormap(colorcet('C06'));
axis image; axis off;
set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
print(h, 'mean_phase_flow_all_sessions2', '-dpdf', '-bestfit', '-painters');