addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%%
load('AP_weights_8points_allsessions.mat');
load('axon_intensity_all_ap.mat');
load('axon_intensity_all_injection.mat');
for i = 1:8
    intensity_all1(:,:,i) = imresize(squeeze(intensity_all(:,:,i)),size(TheColorImage_all,[1,2]));
    injection_intensity1(:,:,i) = imresize(squeeze(injection_intensity(:,:,i)),size(TheColorImage_all,[1,2]));
end

%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
% only select cortex in the atlas
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
%%
point(1,:) = [84,117]; %% SSp-bfd
point(2,:) = [69,122]; %% SSp-n
point(3,:) = [55,118]; %% SSp-m
point(4,:) = [65, 105]; %% SSp-ll
point(5,:) = [73,96]; %% SSp-ul
point(6,:) = [83,93]; %% SSp-tr
point(7,:) = [97 79]; %% RSP
point(8,:) = [115 105]; %% VISp
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% right SSp and MO index
scale = 8;
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
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'right';
[indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; BW_MO(indexMO) =1;
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
%%
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
color1 = colorcet('C06', 'N', 9);
gcamp_mean = mean(TheColorImage_all,4);
figHand = figure('Renderer', 'painters', 'Position', [100 100 600 400]);

axx4 = subplot(1,4,1)
im2 = imagesc(template2);
% im2 = imagesc(not(logical(template2)));
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


ax1 = subplot(1,4,2);
scale = 8;
im2 = imagesc(template2);
% im2 = imagesc(not(logical(template2)));
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
for kkk = 1:8
    TheColorImage = squeeze(gcamp_mean(:,:,kkk));
    % overlay kernel with colormap    
    maxI = max(TheColorImage(:));
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage/maxI;
    colorIntensity = colorIntensity.^2;
    thisColor = color1(kkk,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on
    hb2(kkk) = imshow(thisColorImage);  
    set(hb2(kkk),'AlphaData',colorIntensity);
    axis image; 
end
hold on;
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off
%
ax2 = subplot(1,4,3);
% hb1 = imshow(grayImage,[]);
hb1 = imagesc(template2);
% hb1 = imagesc(not(logical(template2)));
colormap(gray)
set(hb1, 'AlphaData', BW_SSp, 'AlphaDataMapping', 'scaled');
for k = 1:8
    TheColorImage = squeeze(injection_intensity1(:,:,k));
    % overlay kernel with colormap    
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage;
    thisColor = color1(k,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on
    hb2(k) = imshow(thisColorImage);  
    set(hb2(k),'AlphaData',colorIntensity);
    axis image; 
end
set(ax2,'YDir','reverse')
originalSize1 = get(ax2, 'Position');
axis image; 
axis off;
hold on;
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

ax1 = subplot(1,4,4);
scale = 8;
im2 = imagesc(template2);
% im2 = imagesc(not(logical(template2)));
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
for kkk = 1:8
    TheColorImage = squeeze(intensity_all1(:,:,kkk));
    % overlay kernel with colormap    
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage;
    thisColor = color1(kkk,:);
    thisColorImage = cat(3, thisColor(1)*ones(imageSize),...
        thisColor(2)*ones(imageSize), thisColor(3)*ones(imageSize)); 
    hold on
    hb2(kkk) = imshow(thisColorImage);  
    set(hb2(kkk),'AlphaData',colorIntensity);
    axis image; 
end
hold on;
scale3 = 5/8;
hold on;
lineColor = 'k';
lineColor1 = 'w';
hemi = 'right';
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off
%%
print(figHand, 'AP_mean_activity_axon', '-dpdf', '-bestfit', '-painters');
%%
for i = 1:15
    TheColorImage_all1(:,:,:,i) = zscore(TheColorImage_all(:,:,:,i));
end
%%
for kk = 1:15
    dot_real = 0;
    for kkk = 1:8
        gcamp = squeeze(TheColorImage_all1(:,:,kkk,kk));
        axon = squeeze(intensity_all1(:,:,kkk));
        dot_temp = dot(gcamp(:),axon(:));
        dot_real = dot_real+dot_temp;
    end
    dot_real_all(kk,1) = dot_real;
    %%
    for i = 1:1000
        interation = randperm(8);
        intensity_all2 = intensity_all1(:,:,interation);
        dot_perm = 0;
        for kkk = 1:8
            gcamp = squeeze(TheColorImage_all1(:,:,kkk,kk));
            axon = squeeze(intensity_all2(:,:,kkk));
            dot_temp = dot(gcamp(:),axon(:));
            dot_perm = dot_perm+dot_temp;
        end
        dot_perm_all(kk,i) = dot_perm;
    end
    [h1(kk),p(kk)] = ttest2(dot_real_all(kk,1), dot_perm_all(kk,:));
end
%% normalize between 0 and 1
dot_all = [dot_real_all, dot_perm_all];
dot_all = dot_all-min(dot_all(:));
dot_all = dot_all/max(dot_all(:));
%%
save('AP_dot_index_permutation.mat','dot_all','h1','p');
%%
h = figure('Renderer', 'painters', 'Position', [100 100 300 500]);
perm_index = ones(size(1000,1));
for i = 1:15
    scatter(i*perm_index,dot_all(i,2:end),3,'markerFaceColor',[0.5, 0.5,0.5],'markerEdgeColor',[0.5, 0.5,0.5]);
    hold on;
    scatter(i,dot_all(i,1),3,'r','filled');
    mean_perm = mean(dot_all(i,2:end));
    std_perm = std(dot_all(i,2:end));
    hold on;
    scatter(i,mean_perm,28,'k','_');
    hold on;
    errorbar(i,mean_perm,std_perm,'k');
end
xlabel('Mouse ID');
ylabel('Matching index');
%%
print(h, 'AP_dot_index_permutation', '-dpdf', '-bestfit', '-painters');
%%
clear dot_perm_all2
dot_real = 0;
for kkk = 1:8
    gcamp = squeeze(gcamp_mean(:,:,kkk));
    axon = squeeze(intensity_all1(:,:,kkk));
    dot_temp = dot(gcamp(:),axon(:));
    dot_real = dot_real+dot_temp;
end
dot_real_all2 = dot_real;
%
for i = 1:1000
    interation = randperm(8);
    intensity_all2 = intensity_all1(:,:,interation);
    dot_perm = 0;
    for kkk = 1:8
        gcamp = squeeze(gcamp_mean(:,:,kkk));
        axon = squeeze(intensity_all2(:,:,kkk));
        dot_temp = dot(gcamp(:),axon(:));
        dot_perm = dot_perm+dot_temp;
    end
    dot_perm_all2(i,1) = dot_perm;
end
[h2,p2] = ttest2(dot_real_all2, dot_perm_all2);
%% normalize between 0 and 1
dot_all2 = [dot_real_all2; dot_perm_all2];
dot_all2 = dot_all2-min(dot_all2(:));
dot_all2 = dot_all2/max(dot_all2(:));
%%
h = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
edges = [0:0.05:1];
histogram(dot_all2(2:end),edges)
hold on;
xline(dot_all2(1),'r')
xlim([0,1.2]);
ylim([0,150]);
yticks([0,50,100,150]);
yticklabels({'0','50','100','150'});
%%
mean_indx = mean(dot_all2(2:1001));
std_indx = std(dot_all2(2:1001));
sem_indx  =std_indx./sqrt(1000);
[h2,p2] = ttest2(dot_all2(1), dot_all2(2:1001));
%%
print(h, 'AP_mean_dot_index_permutation', '-dpdf', '-bestfit', '-painters');