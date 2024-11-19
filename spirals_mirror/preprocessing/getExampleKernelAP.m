function getExampleKernelAP(T, data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
%%
point(1,:) = [84,117]; %% SSp-bfd
point(2,:) = [69,122]; %% SSp-n
point(3,:) = [55,118]; %% SSp-m
point(4,:) = [65, 105]; %% SSp-ll
point(5,:) = [73,96]; %% SSp-ul
point(6,:) = [83,93]; %% SSp-tr
point(7,:) = [97 79]; %% RSP
point(8,:) = [115 105]; %% VISp
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll',...
    'SSp-ul','SSp-tr','RSP','VISp'};
%%
for kk = 1:15
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    fname1 = [fname '-AP.mat'];
    load(fullfile(data_folder,'spirals_mirror','regression_ap',fname1),...
        'kernel_full2');
    %% loop through 8 example points
    for kkk = 1:8
        TheColorImage_all(:,:,kkk,kk) = squeeze(...
            kernel_full2(:,:,point(kkk,1),point(kkk,2)));        
    end
end
%%
save(fullfile(save_folder,...
    'AP_weight_8points_allsessions.mat'), 'TheColorImage_all');
