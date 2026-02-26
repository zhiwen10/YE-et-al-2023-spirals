function getAxonMapInjection(data_folder,save_folder)
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
dataFolder = fullfile(data_folder,'spirals_mirror\AxonProjectionVolumn');
nameList = ["SSp_bfd"; "SSp_n"; "SSp_m"; "SSp_ll";...
    "SSp_ul"; "SSp_tr"; "RSP"; "VISp"];
nameList1 = strcat(nameList, '_coronal'); 
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
injection_intensity = [];
for k = 1:8
    S = dir(fullfile(dataFolder,[char(nameList1(k)) '*.nrrd']));
    [data, meta] = nrrdread(fullfile(S.folder,S.name));
    data(not(cortexMask)) = 0;
    TheColorImage = squeeze(sum(data,1));
    colorIntensity = TheColorImage/max(TheColorImage(:));
    injection_intensity(:,:,k) = colorIntensity;
end
save(fullfile(save_folder,'axon_intensity_all_injection.mat'),...
    'injection_intensity');