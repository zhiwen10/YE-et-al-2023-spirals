function getAxonMapHEMI(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
%% create cortex mask
[atlas, metaAVGT] = nrrdread(fullfile(data_folder,'tables',...
    'annotation_50.nrrd'));
ctx = '/997/8/567/688/';
cortexMask = create3dMask(ctx,st,atlas);
%%
dataFolder = fullfile(data_folder,'spirals_mirror','AxonProjectionVolumn');
nameList = ["SSp_bfd"; "SSp_n"; "SSp_m"; "SSp_ll";...
    "SSp_ul"; "SSp_tr"; "RSP"; "VISp"];
nameList1 = strcat(nameList, '_coronal'); 
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
hemi = 'left';
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_left = BW_empty; BW_left(indexleft) =1;
%% 
intensity_all = [];
for k = 1:8
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    TheColorImage(~BW_left) = 0;
    colorIntensity = TheColorImage/max(TheColorImage(:));
    intensity_all(:,:,k) = colorIntensity;
end
save(fullfile(save_folder,'axon_intensity_all_hemi.mat'),'intensity_all');