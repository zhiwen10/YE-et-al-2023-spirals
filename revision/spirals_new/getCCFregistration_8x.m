% function getCCFregistration_8x(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
projectedTemplate1 = projectedTemplate1(1:scale:end,1:scale:end);
projectedTemplate1 = projectedTemplate1/max(projectedTemplate1(:));
%% get original and predicted dV by kernel regression.
scale = 8;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
for kk = 1:4
    %%
    clear spiralsT filteredSpirals unique_spirals unique_spirals_unit indx2
    mn = fnames{kk};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T = T_session(T_session.label == "passive",:);
    %%
    session = 5;
    mn = T.MouseID{session};
    tda = T.date(session);
    en = T.folder(session);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    if ~(exist(fullfile(save_folder,fname), 'file') == 2)
        sessionRoot = fullfile(fullfile(data_folder,'task','task_svd',fname));
        mimg = readNPY(fullfile(sessionRoot, 'meanImage.npy'));
        mimg1 = mimg/max(mimg(:));
        mimg1 = mimg1(1:scale:end,1:scale:end);
        clear cpstruct
        h = cpselect(mimg1,projectedTemplate1);
        while ~exist('cpstruct','var')
            pause(2)
        end
        [movingPoints,fixedPoints] = cpstruct2pairs(cpstruct); 
        tform = fitgeotrans(movingPoints,fixedPoints,'affine');    
        save(fullfile(save_folder,fname),'tform');
%        %%
%         sizeTemplate = size(projectedTemplate1);
%         mimgtransformed = imwarp(mimg1,tform,'OutputView',imref2d(sizeTemplate));
%         figure
%         ax2 = subplot(1,1,1); 
%         im = imagesc(mimgtransformed);
%         set(im, 'AlphaData', logical(mimgtransformed), 'AlphaDataMapping', 'scaled');
%         colormap(gray);
%         hold on; 
%         scale1 = 1;
%         overlayOutlines(coords,scale1,'w');
%         set(gca,'Ydir','reverse');
%         axis off; axis image;
    end
end