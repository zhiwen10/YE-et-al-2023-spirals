githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines')))
addpath(genpath(fullfile(githubdir2, 'widefield')))
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
%%
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\spiralDetectionfftn\full_roi';
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35; % frame sampling rate
params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
params.padding = 2*params.halfpadding;
params.th = 1:36:360; 
params.rs = 10:5:20;
params.gridsize = 10;
params.spiralRange = linspace(-pi,pi,5);
params.gsmooth = 0;
params.epochL = 1000;
params.dThreshold = 15;
params.rsRCheck = 10:10:100;
rate = 1;
% only look around roi with gridsize radius
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);
%%
for kk = 14:15
    clear spiralsT unique_spirals_unit filteredSpirals_round 
    mn = T.MouseID{session_all(kk)};
    tda = T.date(session_all(kk));
    en = T.folder(session_all(kk));
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    serverRoot = expPath(mn, td, en);
    [U,V,t,mimg] = get_wf_svd1(serverRoot);
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(roi_folder,[fname '_roi.mat']));
    %%
    mimg = readNPY(fullfile(serverRoot, 'blue','meanImage.npy'));
    params.xsize = size(mimg,1);
    params.ysize = size(mimg,2);
    xsizePadded = params.xsize+params.padding; ysizePadded = params.ysize+params.padding;
    [xx,yy] = meshgrid(min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
        min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
    params.xx = xx; params.yy = yy;
    %% apply mask, this helps speed up spiral detection later
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    mimg2_bw = logical(mimg2);
    %%
    fname1 = [mn '_' tdb '_' num2str(en) '_roi'];
    load(fullfile(roi_folder,[fname1 '.mat']));
    tf = inROI(roi,xx(:),yy(:));
    params.xxRoi = xx(tf); 
    params.yyRoi =yy(tf);
    rate = 1;
    %%
    params.nt = numel(t);
    frameN1 = round((params.nt-70)/100)*100;
    params.frameRange = 36:params.epochL:(36+frameN1-1);
    U1 = U(1:params.downscale:end,1:params.downscale:end,1:50);
    %%
    pwAll = [];
    for kk = 1:numel(params.frameRange)-1
        %%
        tic
        pwAll1 = []; pwAll2 = []; pwAll3 = []; pwAll4 = []; pwAll5 = [];
        frameStart = params.frameRange(kk);
        frameEnd = frameStart+params.epochL-1;
        frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
        dV1 = dV(1:50,frameTemp);
        t1 = t(frameTemp);    
        %%
        % [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_fftn(U1,dV1,t,params,rate);
        [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U1,dV1,t,params,rate);
        tracePhase1 = tracePhase1(:,:,1+35:end-35); % reduce 2*35 frames after filter data 
        tracePhase = padZeros(tracePhase1,params.halfpadding); % pad tracephase with edge zeros
        %%
        nframe = size(tracePhase,3);
        for frame = 1:nframe
            A = squeeze(tracePhase(:,:,frame));
            [pwAll1] = spiralAlgorithm(A,params);
            frameN = frame+frameStart-1;
            pwAll2 = checkClusterXY(pwAll1,params.dThreshold);
            [pwAll3] = doubleCheckSpiralsAlgorithm(A,pwAll2,params);
            [pwAll4] = spatialRefine(A,pwAll3,params);
            [pwAll5] = spiralRadiusCheck2(A,pwAll4,params);
            if not(isempty(pwAll5))
                pwAll5(:,end+1) = frameN; 
            end
            pwAll = [pwAll;pwAll5];
        end
        fprintf('Frame %g/%g; time elapsed %g seconds \n', [frameStart,frameN1, toc])
    end
    %%
    pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;
    frame_count = numel(t);
    save([fname '.mat'], 'pwAll','frame_count');
end