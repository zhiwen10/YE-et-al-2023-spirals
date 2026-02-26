function getSpiralsPhaseMap(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
%% find mask and index for brain region
downscale = 8;
BW = logical(projectedAtlas1);
BW1 = BW(1:downscale:end,1:downscale:end);
[row,col] = find(BW);
brain_index = [col,row];
%%
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
%%
for kk = 1:15
    clear spiralsT filteredSpirals1
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %% load SVD
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    %%
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    %%
    Utransformed = imwarp(U,tform,...
        'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed = Utransformed(1:downscale:end,1:downscale:end,:);
    %% trasnform sprials to atlas space
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);                                          % only use spirals with duration >= 2 frames
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));    
    filteredSpirals(:,1:2) = round(spiralsT); 
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % only use spirals detected within brain ROI
    filteredSpirals = filteredSpirals(lia,:);
    %% only use spirals >70um, single hemisphere, and within SSp ROI
    filteredSpirals1 = filteredSpirals;
    filteredSpirals1 = filteredSpirals(filteredSpirals(:,3)>=70,:);
    filteredSpirals1 = filteredSpirals1(filteredSpirals1(:,4)==1,:);
    location_index = (filteredSpirals1(:,1)>=800 &...
        filteredSpirals1(:,1)<=900 & ...
        filteredSpirals1(:,2)>=500 &...
        filteredSpirals1(:,2)<=650);
    filteredSpirals1 = filteredSpirals1(location_index,:);
    %% get phase maps (current + next frame) for each spiral
    spiral_n = size(filteredSpirals1,1);
    spiral_phase_all = zeros(165,143,2,spiral_n);
    for i = 1:spiral_n
        frameStart = filteredSpirals1(i,5);
        frameEnd = frameStart+100; 
        frameTemp = frameStart-35:frameEnd+35;                             % extra 2*35 frames before filter data 
        dV1 = dV(:,frameTemp);
        [~,~,tracePhase1] = spiralPhaseMap4(...
            Utransformed,dV1,t,params,rate);
        tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate);              % reduce 2*35 frames after filter data 
        tracePhase_current = tracePhase1(:,:,1:2);
        spiral_phase_all(:,:,:,i) = tracePhase_current;        
    end
    %% normalize all sprial phase maps to the same value at point (70,95) 
    spiral_phase_all_norm = zeros(165,143,2,spiral_n);
    for i = 1:spiral_n
        spiral_phase_temp = squeeze(spiral_phase_all(:,:,:,i));
        spiral_phase_temp1 = spiral_phase_temp-spiral_phase_temp(70,95);
        spiral_phase_temp2 = angle(exp(1i*(spiral_phase_temp1))); 
        spiral_phase_all_norm(:,:,:,i) = spiral_phase_temp2;
    end
    %%
    save(fullfile(save_folder,[fname '_mean_flow2']),...
        'spiral_phase_all_norm');
end