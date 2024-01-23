%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
data_folder = 'C:\Users\Steinmetz lab\Documents\git\BIL\function_matching\spirals_70pixels';
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
    load(fullfile(data_folder,[fname '_mean_phase_flow.mat']));
    spiral_phase_all = cat(3,spiral_phase_all,spiral_phase_all_norm);
    vxRaw_all1 = cat(3,vxRaw_all1,vxRaw_all);
    vyRaw_all1 = cat(3,vyRaw_all1,vyRaw_all);
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
[spiral_phase_mean] = circ_mean(spiral_phase_all, [], 3);
%%
skip = 3; zoom_scale = 10;
vxRaw_mean = mean(vxRaw_all1,3); vyRaw_mean= mean(vyRaw_all1,3);
vxRaw2 = nan(size(BW1));
vyRaw2 = nan(size(BW1));
vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
scale2 = 8;
h = figure;
subplot(1,1,1);
im = imagesc(spiral_phase_mean);
colormap(colorcet('C06'));
axis image; axis off;
set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
print(h, 'mean_phase_flow_all_sessions', '-dpdf', '-bestfit', '-painters');