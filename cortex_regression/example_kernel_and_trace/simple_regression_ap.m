%%
% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath(fullfile(githubdir2, 'allenCCF')));
addpath(genpath(fullfile(githubdir2, 'AP_scripts_cortexlab-master')));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
mn = 'ZYE_0012'; td = '2020-10-16'; en = 1;enExp = 2;en = 5;
%%
serverRoot = expPath(mn, td, en);
wf_svd
%%
tda = datetime(td,'InputFormat','yyyy-MM-dd');
formatOut = 'yyyymmdd';
tdb = datestr(tda,formatOut);
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
dfolder = fullfile(folder,mn,td,num2str(en));
fname = [mn '_' tdb '_' num2str(en)];
tformName = dir(fullfile(dfolder,'*tform.mat')).name;
load(fullfile(dfolder,tformName));
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
U = U(:,:,1:50);
V = V(1:50,:);
%%
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
%% only select cortex in the atlas
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'ara_tools-master')))
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
mimg1 = mimg/max(mimg(:));
%%
scale = 8;
%% right SSp index
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[Unew_SSp,Vnew_SSp] = redoSVD(UselectedSSp,V(:,1:58000));
USSp = Unew_SSp(:,1:50);
VSSp = Vnew_SSp(1:50,:);
%% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'right';
[indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[Unew_MO,Vnew_MO] = redoSVD(UselectedMO,V(:,1:58000));
UMO = Unew_MO(:,1:50);
VMO = Vnew_MO(1:50,:);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; BW_MO(indexMO) =1;
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
%% all areas
UtDown = Utransformed(1:scale:end,1:scale:end,1:50);
Utr = reshape(UtDown,size(UtDown,1)*size(UtDown,2),size(UtDown,3));
%% prepare regressor and signal for regression
regressor1 = zscore(VMO,[],2);
regressor = regressor1(:,1:58000);
Vregressor_GPU =gpuArray(regressor);

signal1 = USSp*VSSp(:,1:58000);
signal_GPU = gpuArray(signal1);
%% regression
kk1 = gather(Vregressor_GPU'\signal_GPU');    
k1_real = UMO*kk1;
%% prepare test regressor 
data = UselectedMO*V(:,78000:end);    
% subtract mean as we did before
data = data-mean(data,2);
regressor_test = Unew_MO' * data;   
regressor_test_z = zscore(regressor_test,[],2);
%%
projectedAtlas3 = projectedAtlas1(1:scale:end,1:scale:end);
%% prediction
predicted_signals = kk1'*regressor_test_z(1:50,:);
predicted_signal = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),size(predicted_signals,2));
predicted_signal(indexSSp,:) = predicted_signals;
rawAll = Utr*V(:,78000:end);
BW  = logical(projectedAtlas1(1:scale:end,1:scale:end));
explained_var4 = sseExplainedCal(rawAll ,predicted_signal);
explained_var4 = reshape(explained_var4,size(UtDown,1),size(UtDown,2));
explained_var4(~BW) = 0;
figure; 
im1 = imagesc(explained_var4);
set(im1, 'AlphaData',BW, 'AlphaDataMapping', 'scaled');
axis image; axis off;
colorbar;
caxis([0,1]);
%% kernel project back to pixel space
kernel_full = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),size(projectedAtlas3,1)*size(projectedAtlas3,2));
for j = 1:numel(indexSSp)
     kernel_temp2 = zeros(size(projectedAtlas3,1)*size(projectedAtlas3,2),1);
     kernel_temp2(indexMO) = k1_real(:,j);
     kernel_full(:,indexSSp(j)) = kernel_temp2;
end
kernel_full2 = reshape(kernel_full,[size(projectedAtlas3,1),size(projectedAtlas3,2),...
    size(projectedAtlas3,1),size(projectedAtlas3,2)]);
%%
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
point(8,:) = [110 104]; %% VISp
point(7,:) = [97 79]; %% RSP
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% 3 subareas
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
ctx = '/997/8/567/688/';
%%
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
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; BW_MO(indexMO) =1;
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
hemi = 'right';
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
% mask and Kernel regression map for left
hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);

BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
BW_right = BW_empty; BW_right(indexright) =1;
%%
template1 =zeros(1320,1140);
for k1 = 1:1320
    for k2 = 1:1140
        secondIndex = find(av(k1,:,k2)>1,60,'first');
        if not(isempty(secondIndex))
            template1(k1,k2) = tv(k1,secondIndex(end),k2);
        end
    end
end
template2 = template1(1:8:end,1:8:end);
%% color overlay a cycle
color1 = colorcet('C06', 'N', 9);
mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
h1 = figure;

axx4 = subplot(1,3,1)
im2 = imagesc(mimgtransformed2);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');

hold on;
for kkk = 1:8
    TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
    % overlay kernel with colormap    
    maxI = max(TheColorImage(:));
    % TheGrayscaleImage = squeeze(mimg);
    % hb1 = imshow(TheGrayscaleImage,[]);
    % colormap(ax1,gray);
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage/maxI;
    colorIntensity = colorIntensity.^2;
%     maxIntensity = max(colorIntensity(:));
%     colorIntensity(colorIntensity<0.6*maxIntensity) = 0;
    thisColor = color1(kkk,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on
    hb2(kkk) = imshow(thisColorImage);
    % hold off    
    set(hb2(kkk),'AlphaData',colorIntensity);
    axis image; 

end
hold on;
% overlayOutlines(coords,scale2,'w');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

axx4 = subplot(1,3,2)
im2 = imagesc(mimgtransformed2);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_SSp, 'AlphaDataMapping', 'scaled');
hold on;
for kkk = 1:8 
    scatter(axx4,point(kkk,2),point(kkk,1),12,...
        'MarkerFaceColor',color1(kkk,:),...
        'MarkerEdgeColor','None');
    hold on;
end
scale3 = 5/8;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image;
axis off;

ax2 = subplot(1,3,3);
colormap(ax2,color1(1:8,:));
cb = colorbar('YTickLabel', nameList,'YTick',0.06:0.13:1.05);
cb.TickLabelInterpreter = 'none';
axis off;
originalSize1 = get(ax2, 'Position');
set(ax2, 'Position', [originalSize1(1),originalSize1(1)-0.2,originalSize1(3)*0.6, originalSize1(4)*0.3]);
%%
print(h1, 'KernelMap-ap', '-dpdf', '-bestfit', '-painters');
%%
% pixel(1,:) = [914,933]; %V1
% pixel(2,:) = [590,867]; %S1
% pixel(3,:) = [779,658]; %RSP

pixel(1,:) = [115,105]*8; %V1
pixel(2,:) = [78,118]*8; %S1
pixel(3,:) = [89,78]*8; %RSP

% zye12 anterior cortex conterparts
pixel_a(1,:) = [64,80]*8;
pixel_a(2,:) = [54,88]*8;
pixel_a(3,:) = [72,77]*8;

% zye12 left cortex conterparts
pixel_l(1,:) = [107,41]*8;
pixel_l(2,:) = [75,30]*8;
pixel_l(3,:) = [92,61]*8;

% pixel(1,:) = [111,105]*8; %V1
% pixel(2,:) = [78,118]*8; %S1
% pixel(3,:) = [89,78]*8; %RSP
% ZYE67 anterior cortex conterparts
% pixel(4,:) = [73,84]*8;
% pixel(5,:) = [52,89]*8;
% pixel(6,:) = [78,60]*8;
%% color overlay a cycle
color1 = colorcet('C06', 'N', 9);
mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 400]);
axx4 = subplot(2,3,1)
% im2 = imagesc(mimgtransformed2);
im2 = imagesc(template2);
colormap(gray)
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 8
TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
% overlay kernel with colormap    
maxI = max(TheColorImage(:));
[row,col] = find(TheColorImage==maxI);
% TheGrayscaleImage = squeeze(mimg);
% hb1 = imshow(TheGrayscaleImage,[]);
% colormap(ax1,gray);
imageSize = size(TheColorImage);  
colorIntensity = TheColorImage/maxI;
colorIntensity = colorIntensity.^2;
%     maxIntensity = max(colorIntensity(:));
%     colorIntensity(colorIntensity<0.6*maxIntensity) = 0;
thisColor = color1(kkk,:);
thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
    thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
hold on
hb2(kkk) = imshow(thisColorImage);
hold off

set(hb2(kkk),'AlphaData',colorIntensity);
axis image; 
% hold on;
% scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
hold on;
scatter(axx4,col,row,12,'MarkerFaceColor','k','MarkerEdgeColor','None');
hold on;
% overlayOutlines(coords,scale2,'w');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

axx4 = subplot(2,3,2)
im2 = imagesc(template2);
colormap(gray)
set(im2, 'AlphaData', BW_SSp, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 8;
axis image; 
hold on;
scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

Vraw = V(1:50,78000:end);
traceSSp_predict2 = reshape(predicted_signal,165,143,[]);
subplot(2,3,3);
hold off;
i = 1;
trace_raw1 = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %sensory
trace_raw3 = squeeze(Utransformed(row*scale,col*scale,:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %left hemisphere
trace_predict = squeeze(traceSSp_predict2(round(pixel(i,1)/8),round(pixel(i,2)/8),:))./mimgtransformed(pixel(i,1),pixel(i,2));  % preidicted
plot(1/35:1/35:size(trace_predict,1)/35,trace_predict+0.1,'color',color1(kkk,:));
hold on;
plot(1/35:1/35:size(trace_raw1,2)/35,trace_raw1,'color',[0 0 0 ]);
hold on;
plot(1/35:1/35:size(trace_raw3,2)/35,trace_raw3+0.2,'color',color1(kkk,:));
hold on;
xlim([50,60])
xlabel('Time (s)')
xticks([50:2:60])
xticklabels(string(num2cell(0:2:10)))
ylim([-0.05,0.25]);

% color overlay a cycle
axx4 = subplot(2,3,4);
im2 = imagesc(template2);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 1;
TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2)));
% overlay kernel with colormap    
maxI = max(TheColorImage(:));
[row,col] = find(TheColorImage==maxI);
% TheGrayscaleImage = squeeze(mimg);
% hb1 = imshow(TheGrayscaleImage,[]);
% colormap(ax1,gray);
imageSize = size(TheColorImage);  
colorIntensity = TheColorImage/maxI;
colorIntensity = colorIntensity.^2;
%     maxIntensity = max(colorIntensity(:));
%     colorIntensity(colorIntensity<0.6*maxIntensity) = 0;
thisColor = color1(kkk,:);
thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
    thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
hold on
hb2(kkk) = imshow(thisColorImage);
hold off

set(hb2(kkk),'AlphaData',colorIntensity);
axis image; 
% hold on;
% scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
hold on;
scatter(axx4,col,row,12,'MarkerFaceColor','k','MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

axx4 = subplot(2,3,5)
im2 = imagesc(template2);
colormap(gray)
set(im2, 'AlphaData', BW_right, 'AlphaDataMapping', 'scaled');
hold on;
kkk = 1;
axis image; 
hold on;
scatter(axx4,point(kkk,2),point(kkk,1),12,'MarkerFaceColor',color1(kkk,:),'MarkerEdgeColor','None');
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

Vraw = V(1:50,78000:end);
traceSSp_predict2 = reshape(predicted_signal,165,143,[]);
subplot(2,3,6);
hold off;
i = 2;
trace_raw1 = squeeze(Utransformed(pixel(i,1),pixel(i,2),:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %sensory
trace_raw3 = squeeze(Utransformed(row*scale,col*scale,:))'* Vraw./mimgtransformed(pixel(i,1),pixel(i,2)); %left hemisphere
trace_predict = squeeze(traceSSp_predict2(round(pixel(i,1)/8),round(pixel(i,2)/8),:))./mimgtransformed(pixel(i,1),pixel(i,2));  % preidicted


plot(1/35:1/35:size(trace_predict,1)/35,trace_predict+0.1,'color',color1(kkk,:));
hold on;
plot(1/35:1/35:size(trace_raw1,2)/35,trace_raw1,'color',[0 0 0 ]);
hold on;
plot(1/35:1/35:size(trace_raw3,2)/35,trace_raw3+0.2,'color',color1(kkk,:));

xlim([50,60])
xlabel('Time (s)')
xticks([50:2:60])
xticklabels(string(num2cell(0:2:10)))
ylim([-0.05, 0.25]);
%%
print(h1, 'example_traces_MO2SSp-2', '-dpdf', '-bestfit', '-painters');