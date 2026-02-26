function getSpiralDetectionPredict(T,data_folder, save_folder)
% iterate through all sessions in T for sprial detection
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
    % dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V 
    %% set params for detection
    freq = [2,8];                                                          % data filtering frequency range
    rate = 1;                                                              % set to 1, if no upsampling in time
    params = setSpiralDetectionParams(U,t);                                % set detection parameters
    %% draw brain mask roi, or import from roi folder and detect sprials                                           
    %  main search algorithm 
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'revision2','V_predict',[fname '_v_predict.mat']));
    dV = [zeros(size(V_predict,1),1) diff(V_predict,[],2)]; 
    fname1 = [mn '_' tdb '_' num2str(en) '_roi'];
    load(fullfile(data_folder,'spirals','full_roi',[fname1 '.mat']));
    tf = inROI(roi,params.xx(:),params.yy(:));
    params.xxRoi = params.xx(tf);                                          % only use the grids that inside the roi to save time
    params.yyRoi = params.yy(tf);                                          % only use the grids that inside the roi to save time
    U1 = U(1:params.downscale:end,1:params.downscale:end,1:50);            % no spatial downsampling 
    [pwAll] = spiralDetectionAlgorithm(...                                 % main spiral detection algorithm
        U1,dV,t,params,freq,rate);          
    fname = [mn '_' tdb '_' num2str(en) '_spirals_all.csv'];               % save data file name
    T1 = array2table(pwAll,'VariableNames',...                             % make the sprial matrix into a table with varibale names
        {'spiral_center_x','spiral_center_y',...
        'spiral_radius','spiral_direction','spiral_frame'});            
    writetable(T1,fullfile(save_folder,fname));                            % write to csv file
end