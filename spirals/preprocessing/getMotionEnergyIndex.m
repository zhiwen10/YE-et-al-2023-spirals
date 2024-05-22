function getMotionEnergyIndex(data_folder,save_folder)
%%
% bin spirarity and sync index based on motion energy bins,
% and save data as a .mat file
% bin edges were pre determined based on motion energy histograms across
% sessions

% 'edges': motion energy bin edges
% 'sync','spirality': indexes for each frame in each session
% 'index_sort','sync_sort','spirality_sort': mean index based on motion
% energy bins
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
%%
scale = 1;
count1 = 1;
for kk = [1:6,9:15]                                                        % no motion energy for 7 and 8?
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals\svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root); 
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %% load data
    load(fullfile(data_folder,'spirals\rf_tform',[fname '_tform.mat']));   % load atlas transformation matrix tform;
    load(fullfile(data_folder,'spirals\spirals_index',[fname '_amp.mat']));
    load(fullfile(data_folder,'spirals\spirals_index',...
        [fname '_motion_energy.mat']));
    %%
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    %%
    params.downscale = 8;
    params.lowpass = 0;
    params.gsmooth = 0;
    rate = 1;
    %%
    Utransformed = Utransformed(1:params.downscale:end,...
        1:params.downscale:end,:);
    mimgtransformed = mimgtransformed(1:params.downscale:end,...
        1:params.downscale:end);
    BW = logical(projectedAtlas1);
    BW = BW(1:params.downscale:end,1:params.downscale:end);
    %%
    tStart = 900; tEnd = 1100;                                             % find spirals between time tStart:tEnd
    frameStart = find(t>tStart,1,'first'); 
    frameEnd = find(t>tEnd,1,'first');
    frameTemp = frameStart-35:frameEnd+35;                                 % extra 2*35 frames before filter data 
    dV1 = dV(:,frameTemp);
    [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(...
        Utransformed,dV1,t,params,rate);
    trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
    tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate);                  % reduce 2*35 frames after filter data 
    traceAmp1 = traceAmp1(:,:,1+35/rate:end-35/rate)./mimgtransformed;     % reduce 2*35 frames after filter data 
    %%
    cols = size(tracePhase1,1);
    rows= size(tracePhase1,2);
    %%
    xs = ceil(linspace(-rows/4,rows/4,ceil(rows/2)));
    ys = ceil(linspace(-cols/2,cols/2,cols));
    [Xs,Ys] = meshgrid(xs,ys);
    anglein = atan2(Ys,Xs);                                                % expected angle for centered sprial
    BW_half = BW(:,72:end);
    %%
    for frame = 1:size(tracePhase1,3)
        phase_hemi = squeeze(tracePhase1(:,72:end,frame));
        phase_hemi(not(BW_half)) = nan;
        pixel_count = sum(BW_half(:));
        sync(frame,count1)  = abs(sum(exp(j*phase_hemi(:)),"omitnan"))...
            ./pixel_count;
        phase_diff = phase_hemi-anglein;
        spirality(frame,count1) = abs(sum(exp(j*phase_diff(:)),"omitnan"))...
            ./pixel_count;
        amp_hemi = squeeze(traceAmp1(:,72:end,frame));
        amp_hemi(not(BW_half)) = nan;
        amp_hemi = amp_hemi(:);
        amp_mean(frame,count1) = mean(amp_hemi,"omitnan");
    end
    current_energy = image_energy2(frameStart:frameEnd);
    energy_all(:,count1) = current_energy;

    count1 = count1+1;
end
%%
for kk = 1:13
    sync1 = sync(:,kk);
    spirality1 = spirality(:,kk);
    energy1 = energy_all(:,kk);
    index_unity1 = sqrt(sync1.^2+spirality1.^2);
    index_unity(:,kk) = index_unity1;
    edges = [0:50000:1000000];
    clear indx1
    for i = 1:numel(edges)-1
        clear a
        a = find(energy1>=edges(i) & energy1<edges(i+1));
        indx1{i} = a;
        if numel(a)>10
            sync_sort(i,kk) = mean(sync1(indx1{i}));
            spirality_sort(i,kk) = mean(spirality1(indx1{i}));
            index_sort(i,kk) = mean(index_unity1(indx1{i}));
        else
            sync_sort(i,kk) = nan;
            spirality_sort(i,kk) = nan;
            index_sort(i,kk) = nan;
        end
    end
end
%%
save(fullfile(save_folder,'energy_sync_bin_0_0.01_all.mat'),...
    'edges','sync','spirality','index_sort','sync_sort','spirality_sort');
