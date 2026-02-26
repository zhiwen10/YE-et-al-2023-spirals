function h3h = plotKernelMapsHEMI(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
load(fullfile(data_folder,'spirals_mirror','matching_index',...
    'hemi_weights_8points_allsessions.mat'));
load(fullfile(data_folder,'spirals_mirror','matching_index',...
    'axon_intensity_all_hemi.mat'));
load(fullfile(data_folder,'spirals_mirror','matching_index',...
    'axon_intensity_all_injection.mat'));
for i = 1:8
    intensity_all1(:,:,i) = imresize(squeeze(...
        intensity_all(:,:,i)),size(TheColorImage_all,[1,2]));
    injection_intensity1(:,:,i) = imresize(squeeze(...
        injection_intensity(:,:,i)),size(TheColorImage_all,[1,2]));
end
%%
spath = string(st.structure_id_path);
point(1,:) = [84,117]; %% SSp-bfd
point(2,:) = [69,122]; %% SSp-n
point(3,:) = [55,118]; %% SSp-m
point(4,:) = [65, 105]; %% SSp-ll
point(5,:) = [73,96]; %% SSp-ul
point(6,:) = [83,93]; %% SSp-tr
point(7,:) = [97 79]; %% RSP
point(8,:) = [115 105]; %% VISp
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ul','SSp-ll','SSp-tr','RSP','VISp'};
%%
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
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);

hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);

BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
BW_right = BW_empty; BW_right(indexright) =1;
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
%% mask and Kernel regression map for right and left
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
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
% mask and Kernel regression map for left
hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);

BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
BW_right = BW_empty; BW_right(indexright) =1;
%% right SSp and MO index
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
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'right';
[indexMO,UselectedMO] = select_area(frotalArea,spath,st,...
    coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; BW_MO(indexMO) =1;
BW_SSp = BW_empty; BW_SSp(indexSSp) =1;
%%
color1 = colorcet('C06', 'N', 9);
gcamp_mean = mean(TheColorImage_all,4);
h3h = figure('Renderer', 'painters', 'Position', [100 100 500 400]);
axx4 = subplot(1,4,1);
im2 = imagesc(template1);
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
im2 = imagesc(template1);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_left, 'AlphaDataMapping', 'scaled');
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
hemi = 'left';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off
%
ax2 = subplot(1,4,3);
hb1 = imagesc(template1);
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
im2 = imagesc(template1);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_left, 'AlphaDataMapping', 'scaled');
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
hemi = 'left';
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off
%%
print(h3h, fullfile(save_folder,'Fig3h_hemi_mean_activity_axon'),...
    '-dpdf', '-bestfit', '-painters');