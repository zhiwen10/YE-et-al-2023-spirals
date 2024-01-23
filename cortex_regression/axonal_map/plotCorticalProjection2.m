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
nameList = ["SSp_bfd"; "SSp_n"; "SSp_m"; "SSp_ll"; "SSp_ul"; "SSp_tr"; "RSP"; "VISp"];
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
projectedAtlas2 = projectedAtlas1(1:5:end,1:5:end);
BW = logical(projectedAtlas2);
%% plot
color1 = colorcet('C06', 'N', 9);
figHand = figure;
ax1 = subplot(1,2,1);
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

hb1 = imshow(grayImage,[]);
for k = 1:8
    clear data row col
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    imageSize = size(TheColorImage);
    colorIntensity = TheColorImage/max(TheColorImage(:));
    thisColor = color1(k,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on;
    hb2(k) = imshow(thisColorImage);
    hold off;
    set(hb2(k),'AlphaData',colorIntensity);
end
set(ax1,'YDir','reverse')
originalSize1 = get(ax1, 'Position');
axis image; 
axis off;
hold on;
% plotOutline(root1,st,atlas1,[]);
% plotOutline(ctx,st,atlas1,[]);
% plotOutline(ctx,st,atlas1,-1);
scale3 = 1;
plotOutline(maskPath(1:3),st,atlas1,[],scale3);
plotOutline(maskPath(4),st,atlas1,[],scale3);
plotOutline(maskPath(5),st,atlas1,[],scale3);
plotOutline(maskPath(6:11),st,atlas1,[],scale3);
text(80, -50, 'Projection Map','FontSize',12); 
set(hb1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
%
ax2 = subplot(1,2,2);
colormap(ax2,color1);
cb = colorbar('YTickLabel', nameList,'YTick',0.06:0.13:1.05);
cb.TickLabelInterpreter = 'none';
axis off;
set(ax2, 'Position', [originalSize1(1)+0.2,originalSize1(1)+0.14,originalSize1(3)*0.6, originalSize1(4)*0.6]);
%
print(figHand, 'projection_map', '-dpdf', '-bestfit', '-painters');