gcamp_mean = mean(TheColorImage_all,4);
%%
load('axon_intensity_all_ap.mat');
for i = 1:8
    intensity_all1(:,:,i) = imresize(squeeze(intensity_all(:,:,i)),size(TheColorImage_all,[1,2]));
end
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
figure;
ax1 = subplot(1,2,1);
scale = 8;
mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
im2 = imagesc(mimgtransformed2);
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
%%
ax1 = subplot(1,2,2);
scale = 8;
mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
im2 = imagesc(mimgtransformed2);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');
for kkk = 1:8
    TheColorImage = squeeze(intensity_all1(:,:,kkk));
    % overlay kernel with colormap    
    maxI = max(TheColorImage(:));
    imageSize = size(TheColorImage);  
    colorIntensity = TheColorImage/maxI;
    % colorIntensity = colorIntensity.^2;
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
dot_real = 0;
for kkk = 1:8
    gcamp = squeeze(gcamp_mean(:,:,kkk));
    axon = squeeze(intensity_all1(:,:,kkk));
    dot_temp = dot(gcamp(:),axon(:));
    dot_real = dot_real+dot_temp;
end
%%
for i = 1:100
    interation = randperm(8);
    intensity_all2 = intensity_all1(:,:,interation);
    dot_perm = 0;
    for kkk = 1:8
        gcamp = squeeze(gcamp_mean(:,:,kkk));
        axon = squeeze(intensity_all2(:,:,kkk));
        dot_temp = dot(gcamp(:),axon(:));
        dot_perm = dot_perm+dot_temp;
    end
    dot_perm_all(i,1) = dot_perm;
end