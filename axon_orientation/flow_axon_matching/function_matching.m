%%
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\nrrdread'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
kk = 6;
mn = T.MouseID{session_all(kk)};
tda = T.date(session_all(kk));
en = T.folder(session_all(kk));
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
serverRoot = expPath(mn, td, en);
%% SVD, plot example trace, overlay with pupil and stim time
serverRoot = expPath(mn, td, en);
[U,V,t,mimg] = get_wf_svd1(serverRoot);
dV = [zeros(size(V,1),1) diff(V,[],2)];
%% get atlas mask and outlines
% load coords for atlas
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
% load projectedAtlas and projectedTemplate
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
dfolder = fullfile(folder,mn,td,num2str(en));
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(dfolder,[fname '_spiral_group.mat']));
tformName = dir(fullfile(dfolder,'*tform.mat')).name;
load(fullfile(dfolder,tformName));
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
BW = logical(projectedAtlas1);
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
%%
[row,col] = find(BW);
brain_index = [col,row];
%%
clear spiralsT
spiral_duration = cellfun(@(x) size(x,1), archiveCell);
indx2 = (spiral_duration>=3);
groupedCells = archiveCell(indx2);
filteredSpirals = cell2mat(groupedCells);   
[spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,filteredSpirals(:,1),filteredSpirals(:,2));    
filteredSpirals(:,1:2) = round(spiralsT); 
[lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');
filteredSpirals = filteredSpirals(lia,:);
%%
clear filteredSpirals1
filteredSpirals1 = filteredSpirals(filteredSpirals(:,3)>=70,:);
filteredSpirals1 = filteredSpirals1(filteredSpirals1(:,4)==1,:);
location_index = (filteredSpirals1(:,1)>=800&filteredSpirals1(:,1)<=900&filteredSpirals1(:,2)>=500&filteredSpirals1(:,2)<=650);
filteredSpirals1 = filteredSpirals1(location_index,:);
%%
frame_count = max(filteredSpirals(:,5));
hist_bin = 40;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
figure;
ax1 = subplot(1,1,1);
[unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(filteredSpirals1,hist_bin);
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_count*35; % spirals/(mm^2*s)
scatter(unique_spirals(:,1),unique_spirals(:,2),1,unique_spirals_unit,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
colormap(ax1,hot);
axis off; axis image;
xlim(ax1,[0,1140]);
ylim(ax1,[0,1320]);
%% sanity check
useGPU = 0;
spiral_n = size(filteredSpirals1,1);
for i = 10:20
    frameStart = filteredSpirals1(i,5);
    frameEnd = frameStart+100; 
    frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
    dV1 = dV(:,frameTemp);
    [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(Utransformed,dV1,t,params,rate);
    trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
    tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
    tracePhase_current = tracePhase1(:,:,1:2);
    tracePhase_current = permute(tracePhase_current,[3 1 2]);
    [vxRaw,vyRaw] = HS_flowfield(tracePhase_current,useGPU);
    figure; 
    ax1 = subplot(1,1,1);
    im = imagesc(squeeze(tracePhase1(:,:,1)));
    colormap(ax1,colorcet('C06'));
    axis image; axis off;
    set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
    hold on;
    scatter(filteredSpirals1(i,1)/8,filteredSpirals1(i,2)/8,16,'k','filled');
    hold on;
    vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
    vxRaw2 = nan(size(vxRaw));
    vyRaw2 = nan(size(vyRaw));
    skip = 3; zoom_scale = 2;
    vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
    hold on;
    imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
end
%%
useGPU = 0;
spiral_n = size(filteredSpirals1,1);
vxRaw_all = []; vyRaw_all =[];
vxRaw_all = zeros(165,143,spiral_n); vyRaw_all= zeros(165,143,spiral_n);
spiral_phase_all = zeros(165,143,spiral_n);
for i = 1:spiral_n
    tracePhase_current  = [];
    vxRaw = [];vyRaw = [];
    frameStart = filteredSpirals1(i,5);
    frameEnd = frameStart+100; 
    frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
    dV1 = dV(:,frameTemp);
    [~,~,tracePhase1] = spiralPhaseMap4(Utransformed,dV1,t,params,rate);
    tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
    tracePhase_current = tracePhase1(:,:,1:2);
    spiral_phase_all(:,:,i) = squeeze(tracePhase_current(:,:,1));
    tracePhase_current = permute(tracePhase_current,[3 1 2]);    
    [vxRaw,vyRaw] = HS_flowfield(tracePhase_current,useGPU);
    vxRaw = squeeze(vxRaw);vyRaw = squeeze(vyRaw);
    vxRaw_all(:,:,i) = vxRaw;
    vyRaw_all(:,:,i) = vyRaw;   
end
%%
vxRaw_mean = mean(vxRaw_all,3); vyRaw_mean= mean(vyRaw_all,3);
%%
figure;
ax1 = subplot(1,1,1);
% im_phase = imagesc(framea);
% colormap(ax1,colorcet('C06'));
% hold on;
vxRaw_mean = squeeze(vxRaw_mean);
vyRaw_mean = squeeze(vyRaw_mean);
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
skip = 3; zoom_scale = 10;
vxRaw2(1:skip:end,1:skip:end) = vxRaw_mean(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw_mean(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
% set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
scale2 = 8;
overlayOutlines(coords,scale2);
axis image; axis off;
set(gca,'Ydir','reverse')
%%
figure;
scale2 = 8;
overlayOutlines(coords,scale2);
axis image;
set(gca,'Ydir','reverse')
ssp_bfd_point = [70,95];
scatter(ssp_bfd_point(2),ssp_bfd_point(1),'r','filled');
%%
spiral_phase_all_norm = zeros(165,143,spiral_n);
for i = 1:spiral_n
    spiral_phase_temp = squeeze(spiral_phase_all(:,:,i));
    spiral_phase_temp1 = spiral_phase_temp-spiral_phase_temp(70,95);
    spiral_phase_temp2 = angle(exp(1i*(spiral_phase_temp1))); 
    spiral_phase_all_norm(:,:,i) = spiral_phase_temp2;
end
%%
[spiral_phase_mean] = circ_mean(spiral_phase_all_norm, [], 3);
%%
h = figure;
subplot(1,1,1);
im = imagesc(spiral_phase_mean);
colormap(colorcet('C06'));
axis image; axis off;
set(im, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
print(h, 'mean_phase_flow', '-dpdf', '-bestfit', '-painters');
%%
flow_angle = round(rad2deg(atan(vyRaw_mean(:)./vxRaw_mean(:))));
colora2= colorcet('C06','N',180);
x1 = linspace(-90,90,180);
colora3 = interp1(x1,colora2,flow_angle);
colora3 = reshape(colora3,[165,143,3]);
%%
[row,col] = find(not(isnan(vyRaw2)));
figure;
scale2 = 8;
overlayOutlines(coords,scale2);
scale1 = 1;
for i = 1:numel(row)
    plot([col(i)-vxRaw2(row(i),col(i))*scale1,col(i)+vxRaw2(row(i),col(i))*scale1],...
        [row(i)-vyRaw2(row(i),col(i))*scale1,row(i)+vyRaw2(row(i),col(i))*scale1],...
        'color',colora3(row(i),col(i),:),'LineWidth',1);
    hold on;
end
set(gca,'Ydir','reverse')