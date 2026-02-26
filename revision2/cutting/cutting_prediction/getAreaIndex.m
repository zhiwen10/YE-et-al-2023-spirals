function area_index = getAreaIndex(data_folder)
%% load atlas brain horizontal projection and outline
% data_folder = 'E:\spiral_data_share\data'; 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
spath = string(st.structure_id_path);
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
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
hemi = [];
scale = 1;
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
hemi = [];
[indexMO,UselectedMO] = select_area(...
    frotalArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO_L = BW_empty; BW_MO_L(indexMO) = 1;
BW_SSp_L = BW_empty; BW_SSp_L(indexSSp) = 1;
BW_MO_L(:,570:end) = 0; BW_SSp_L(:,570:end) = 0;
[row1,col1] = find(BW_MO_L); MO_L_index = [col1,row1];
[row2,col2] = find(BW_SSp_L); SSp_L_index = [col2,row2];

BW_MO_R = BW_empty; BW_MO_R(indexMO) = 1;
BW_SSp_R = BW_empty; BW_SSp_R(indexSSp) = 1;
BW_MO_R(:,1:570) = 0; BW_SSp_R(:,1:570) = 0;
[row1,col1] = find(BW_MO_R); MO_R_index = [col1,row1];
[row2,col2] = find(BW_SSp_R); SSp_R_index = [col2,row2];

area_index = {MO_L_index,MO_R_index,SSp_L_index,SSp_R_index};