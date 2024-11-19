function getSpiralDetectionFftnRaw(T,freq,data_folder,save_folder)
% iterate through all sessions in T for sprial detection
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
for kk = 1:size(T,1)
    %% session info
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read svd components (U,V,t) from processed data folder
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    %% set params for detection
    rate = 1;                                                              % set to 1, if no upsampling in time
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','fft_roi',[fname '_roi.mat']));                         
    params = setSpiralDetectionParamsFFTN(mimg,t,roi_ap,roi_ml);           % set detection parameters, with square ROI
    U2 = U(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2),:);
    U1 = U2(1:params.downscale:end,1:params.downscale:end,1:50);           % only use U space within square ROI
    %%
    [pwAll] = spiralDetectionAlgorithm(U1,dV,t,params,freq,rate);          % same detection algorithm as the main one
    pwAll1 = pwAll;
    pwAll1(:,1) = roi_ml(1)+pwAll1(:,1)-1;
    pwAll1(:,2) = roi_ap(1)+pwAll1(:,2)-1;
    frame_count = numel(t);
    save(fullfile(save_folder,freq_folder,'control',[fname '.mat']), ...
        'pwAll1','roi_ml','roi_ap','frame_count');
end
