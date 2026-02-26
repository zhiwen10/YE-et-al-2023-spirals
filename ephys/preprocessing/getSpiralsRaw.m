function getSpiralsRaw(T,data_folder,save_folder)
%%
for kk = 1:size(T,1)
    clear pwAll1 params U V t mimg dV roi U1
    ops = get_session_info2(T,kk,data_folder);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    subfolder = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    session_root = fullfile(data_folder,'ephys','svd_spikes',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    %% set params for detection
    freq = [2,8];                                                          % data filtering frequency range
    rate = 1;                                                              % set to 1, if no upsampling in time
    params = setSpiralDetectionParams(U,t);                                % set detection parameters
    %%
    fname1 = [fname '_roi'];
    load(fullfile(data_folder,'ephys','roi',[fname1 '.mat']));               % only detect spirals within ROI
    tf = inROI(roi,params.xx(:),params.yy(:));
    params.xxRoi = params.xx(tf);                                          % only use the grids that inside the roi to save time
    params.yyRoi = params.yy(tf);                                          % only use the grids that inside the roi to save time
    U1 = U(1:params.downscale:end,1:params.downscale:end,1:50);            % no spatial downsampling 
    [pwAll] = spiralDetectionAlgorithm(U1,dV,t,params,freq,rate);          % main spiral detection algorithm
    %%
    fname2 = [fname '_spirals.mat'];
    save(fullfile(save_folder,fname2),'pwAll');
end

