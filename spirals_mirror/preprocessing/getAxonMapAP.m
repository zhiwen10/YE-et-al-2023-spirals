function getAxonMapAP(data_folder,save_folder)
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
%% MO index
Utransformed = zeros(size(projectedAtlas1));
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = 'right';
scale = 5;
[indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO = BW_empty; 
BW_MO(indexMO) =1;
%%
intensity_all = [];
for k = 1:8
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));           % load 3d axon projection volumn 
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    TheColorImage(~BW_MO) = 0;
    colorIntensity = TheColorImage/max(TheColorImage(:));
    intensity_all(:,:,k) = colorIntensity;
end
save(fullfile(save_folder,'axon_intensity_all_ap.mat'),'intensity_all');