function getBrainMask(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%% only select cortex in the atlas
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas1,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
scale = 8;
template2 = projectedTemplate1(1:scale:end,1:scale:end);
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
ctx = '/997/8/567/688/';
% mask and Kernel regression map for right and left
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
%% left and right SSp
hemi = 'right';
[indexright,UselectedRight] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
% mask and Kernel regression map for left
hemi = 'left';
[indexleft,UselectedLeft] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);

BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_SSp_left = BW_empty; BW_SSp_left(indexleft) =1;
BW_SSp_right = BW_empty; BW_SSp_right(indexright) =1;
%% left and right MO
hemi = 'right';
[indexSSp,UselectedSSp] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'left';
[indexMO_left,UselectedMO_left] = select_area(...
    frotalArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);

hemi = 'right';
[indexMO,UselectedMO] = select_area(...
    frotalArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO_left = BW_empty; BW_MO_left(indexMO_left) =1;
BW_MO_right = BW_empty; BW_MO_right(indexMO) =1;
%% plot division
color1 = colorcet('C06', 'N', 9);
h3b = figure('Renderer', 'painters', 'Position', [100 100 700 700]);
ax1 = subplot(2,2,1);
hb1 = imagesc(template2);
colormap(gray)
set(hb1, 'AlphaData', BW_MO_left, 'AlphaDataMapping', 'scaled');
set(ax1,'YDir','reverse')
axis image; axis off;
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'left';
hold on;
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off

ax2 = subplot(2,2,2);
hb1 = imagesc(template2); 
colormap(gray)
set(hb1, 'AlphaData', BW_MO_right, 'AlphaDataMapping', 'scaled');
set(ax1,'YDir','reverse')
axis image; axis off;
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'right';
hold on;
plotOutline(maskPath(1:3),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis on; axis off

% hemi_division
ax3 = subplot(2,2,3);
hb1 = imagesc(template2);
colormap(gray)
set(hb1, 'AlphaData', BW_SSp_left, 'AlphaDataMapping', 'scaled');
set(ax3,'YDir','reverse')
axis image; axis off;
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'left';
hold on;
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis image; axis off;

ax4 = subplot(2,2,4);
hb1 = imagesc(template2);
colormap(gray)
set(hb1, 'AlphaData', BW_SSp_right, 'AlphaDataMapping', 'scaled');
set(ax4,'YDir','reverse')
axis image; axis off;
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = 'right';
hold on;
plotOutline(maskPath(4),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(5),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(6:11),st,atlas1,hemi,scale3,lineColor1);
plotOutline(maskPath(4:11),st,atlas1,hemi,scale3,lineColor);
set(gca,'Ydir','reverse')
axis on; axis off
%%
save(fullfile(save_folder,'area_mask.mat'),'BW_SSp_left','BW_SSp_right',...
    'BW_MO_left','BW_MO_right');

    
