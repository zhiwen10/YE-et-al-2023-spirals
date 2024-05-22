function getROI(T,data_folder, save_folder)
% iterate through all sessions in T for drawing brain roi
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    %% set params for detection
    params = setSpiralDetectionParams(U,t);                                % set detection parameters
    %% draw brain mask roi
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    figure; 
    ax1 = imagesc(mimg2);
    roi = drawpolygon;
    fname1 = [mn '_' tdb '_' num2str(en) '_roi.mat'];
    save(fullfile(save_folder,fname1),'roi');                                            % results saved in <full_roi> folder 
end