function getMapsSession(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
scale = 8;
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
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
%% MO
ffolder = fullfile(data_folder,'spirals_mirror','regression_ap');
colorIntensityAll = zeros(size(BW,1),size(BW,2),8,15);
mimgtransformed2All = zeros(size(BW,1),size(BW,2),15);
count = 1;
for kk = 1:15                                                                    % use ZYE_0012 as an example session
    mn = T.MouseID{kk};
    en = T.folder(kk);    
    tda = T.date(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform
    %% read svd components (U,V,t) from processed data folder
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    mimg = readNPY(fullfile(session_root, 'meanImage.npy')); 
    fname1 = [fname '-AP.mat'];
    load(fullfile(ffolder,fname1));
    %%
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
    mimgtransformed2All(:,:,kk) = mimgtransformed2;
    for kkk = 1:8
        TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2))); 
        maxI = max(TheColorImage(:));
        colorIntensity = TheColorImage/maxI;
        colorIntensity = colorIntensity.^2;
        colorIntensityAll(:,:,kkk,kk) = colorIntensity;
    end
    count = count+1;
end
save(fullfile(save_folder,'kernelMaps_allSessions_AP.mat'),...
    'colorIntensityAll','mimgtransformed2All');
%% LEFT TO RIGHT
clear colorIntensityAll mimgtransformed2All
ffolder = fullfile(data_folder,'spirals_mirror','regression_hemi');
colorIntensityAll = zeros(size(BW,1),size(BW,2),8,15);
mimgtransformed2All = zeros(size(BW,1),size(BW,2),15);
count = 1;
for kk = 1:15                                                                    % use ZYE_0012 as an example session
    mn = T.MouseID{kk};
    en = T.folder(kk);    
    tda = T.date(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform
    %% read svd components (U,V,t) from processed data folder
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    mimg = readNPY(fullfile(session_root, 'meanImage.npy')); 
    fname1 = [fname '-hemi.mat'];
    load(fullfile(ffolder,fname1));
    %%
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    mimgtransformed2 = mimgtransformed(1:scale:end,1:scale:end);
    mimgtransformed2All(:,:,kk) = mimgtransformed2;
    for kkk = 1:8
        TheColorImage = squeeze(kernel_full2(:,:,point(kkk,1),point(kkk,2))); 
        maxI = max(TheColorImage(:));
        imageSize = size(TheColorImage);  
        colorIntensity = TheColorImage/maxI;
        colorIntensity = colorIntensity.^2;
        colorIntensityAll(:,:,kkk,kk) = colorIntensity;
    end
    count = count+1;
end
save(fullfile(save_folder,'kernelMaps_allSessions_hemi.mat'),...
    'colorIntensityAll','mimgtransformed2All');
