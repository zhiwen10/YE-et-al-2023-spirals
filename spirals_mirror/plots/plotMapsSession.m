function [hs11a, hs11b] = plotMapsSession(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%%
point(8,:) = [115 105]; %% VISp
point(7,:) = [97 79]; %% RSP
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% cortex surface outline
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
scale = 8;
%% mask and Kernel regression map for right and left
spath = string(st.structure_id_path);
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
%% MO
load(fullfile(data_folder,'spirals_mirror','regression_kernels',...
    'kernelMaps_allSessions_AP.mat'));
%%
hs11a = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
count = 1;
color1 = colorcet('C06', 'N', 9);
axx4 = subplot(4,4,1);
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

for kk = 1:15                                                                    % use ZYE_0012 as an example session
    %% color overlay a cycle
    mimgtransformed2 = squeeze(mimgtransformed2All(:,:,kk));
    axx4 = subplot(4,4,kk+1);
    im2 = imagesc(mimgtransformed2);
    colormap(gray)
    BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
    set(im2, 'AlphaData', BW_MO, 'AlphaDataMapping', 'scaled');

    hold on;
    for kkk = 1:8
        colorIntensity = squeeze(colorIntensityAll(:,:,kkk,kk));
        imageSize = size(colorIntensity);  
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
    title(T.MouseID{kk}, 'fontsize',6,'Interpreter', 'none');
    count = count+1;
end
print(hs11a, fullfile(save_folder,'FigS11a_map_AP.pdf'),...
    '-dpdf', '-bestfit', '-painters');
%% LEFT TO RIGHT
load(fullfile(data_folder,'spirals_mirror','regression_kernels',...
    'kernelMaps_allSessions_hemi.mat'));
hs11b = figure('Renderer', 'painters', 'Position', [100 100 600 800]);
axx4 = subplot(4,4,1);
im2 = imagesc(template1);
colormap(gray)
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
set(im2, 'AlphaData', BW_right, 'AlphaDataMapping', 'scaled');
hold on;
for kkk = 1:8 
    scatter(axx4,point(kkk,2),point(kkk,1),12,...
        'MarkerFaceColor',color1(kkk,:),...
        'MarkerEdgeColor','None');
    hold on;
end
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'right';
hold on;
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off
count = 1;
color1 = colorcet('C06', 'N', 9);

for kk = 1:15                                                                    % use ZYE_0012 as an example session
    %% color overlay a cycle
    mimgtransformed2 = squeeze(mimgtransformed2All(:,:,kk));
    axx4 = subplot(4,4,kk+1);
    im2 = imagesc(mimgtransformed2);
    colormap(gray)
    BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
    set(im2, 'AlphaData', BW_left, 'AlphaDataMapping', 'scaled');

    hold on;
    for kkk = 1:8
        colorIntensity = squeeze(colorIntensityAll(:,:,kkk,kk));
        imageSize = size(colorIntensity); 
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
    lineColor = 'k'; lineColor1 = 'w';
    hemi = 'left';
    hold on;
    plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
    plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
    plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
    plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
    set(gca,'Ydir','reverse')
    axis image; axis off
    title(T.MouseID{kk}, 'fontsize',6,'Interpreter', 'none');
    count = count+1;
end
print(hs11b, fullfile(save_folder,'FigS11b_map_hemi.pdf'),...
    '-dpdf', '-bestfit', '-painters');
