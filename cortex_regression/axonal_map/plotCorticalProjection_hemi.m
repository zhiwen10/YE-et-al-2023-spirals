%% addpath
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'matGeom\matGeom')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection'));
%% scale down brainGridData
fpath = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions';
brainGridData = readNPY(fullfile(fpath, 'brainGridData.npy'));
brainGridData = brainGridData(:,[3,1,2]);
brainGridData = brainGridData/4.7188;
%% import data
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
nameList = ["SSp_bfd"; "SSp_n"; "SSp_m"; "SSp_ul"; "SSp_ll"; "SSp_tr"; "RSP"; "VISp"];
nameList1 = strcat(nameList, '_coronal');
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
ctx = '/997/8/567/688/';
cortexMask = create3dMask(ctx,st,atlas);
%%
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
[template, meta] = nrrdread(fullfile(dataFolder,'average_template_50.nrrd'));   
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
%%
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
%%
projectedAtlas2 = projectedAtlas1(1:5:end,1:5:end);
BW = logical(projectedAtlas2);
%% only select cortex in the atlas
githubDir1 = 'C:\Users\Steinmetz lab\Documents\MATLAB';
addpath(genpath(fullfile(githubDir1, 'ara_tools-master')))
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
tv = readNPY(fullfile(dfolder, 'template_volume_10um.npy')); % grey-scale "background signal intensity" 
av = readNPY(fullfile(dfolder,'annotation_volume_10um_by_index.npy')); % the number at each pixel labels the area, see note below
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
spath = string(st.structure_id_path);
%% mask and Kernel regression map for right
% sensoryArea spath projectedAtlas1 projectedTemplate1 scale Utransformed V(:,1:58000)
% UselectedRight Unew_right,Vnew_right
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
scale = 5;
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
BW_right = BW_empty; BW_right(indexright) =1;
%%
% plot each brain injection vloumn
grayImage = zeros(264,228);
for i1 = 1:264
    for i2 = 1:228
        temp = squeeze(atlas(:,i1,i2));
        if any(temp)
            index = find(temp>0,10);
            indexN = numel(index);
            grayImage(i1,i2) = template(index(indexN),i1,i2);
        end
    end
end
%% plot
color1 = colorcet('C06', 'N', 9);
figHand = figure;
ax1 = subplot(1,3,1);
% hb1 = imshow(grayImage,[]);
hb1 = imagesc(grayImage);
colormap(gray)
set(hb1, 'AlphaData', BW_left, 'AlphaDataMapping', 'scaled');
intensity_all = [];
for k = 1:8
    clear data row col
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    TheColorImage(~BW_left) = 0;
    imageSize = size(TheColorImage);
    colorIntensity = TheColorImage/max(TheColorImage(:));
    thisColor = color1(k,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on;
    hb2(k) = imshow(thisColorImage);
    hold off;
    set(hb2(k),'AlphaData',colorIntensity);
    intensity_all(:,:,k) = colorIntensity;
end
set(ax1,'YDir','reverse')
originalSize1 = get(ax1, 'Position');
axis image; 
axis off;
hold on;

scale3 = 1;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'left';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image;
axis off;
%
ax2 = subplot(1,3,2);
% hb1 = imshow(grayImage,[]);
hb1 = imagesc(grayImage);
colormap(gray)
set(hb1, 'AlphaData', BW_right, 'AlphaDataMapping', 'scaled');
for k = 1:8
    clear data row col
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    TheColorImage(~BW_right) = 0;
    imageSize = size(TheColorImage);
    colorIntensity = TheColorImage/max(TheColorImage(:));
    thisColor = color1(k,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on;
    hb2(k) = imshow(thisColorImage);
    hold off;
    set(hb2(k),'AlphaData',colorIntensity);
    injection_intensity(:,:,k) = colorIntensity;
end
set(ax2,'YDir','reverse')
originalSize1 = get(ax2, 'Position');
axis image; 
axis off;
hold on;

scale3 = 1;
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

ax3 = subplot(1,3,3);
colormap(ax3,color1(1:8,:));
cb = colorbar('YTickLabel', nameList,'YTick',0.06:0.13:1.05);
cb.TickLabelInterpreter = 'none';
axis off;
set(ax3, 'Position', [originalSize1(1)+0.3,originalSize1(1),originalSize1(3)*0.6, originalSize1(4)*0.3]);
%%
print(figHand, 'projection_map_left_right', '-dpdf', '-bestfit', '-painters');