githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')));
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield')));
addpath(genpath(fullfile(githubdir2, 'npy-matlab')));
%% get original and predicted dV by kernel regression.
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
T = readtable('session_list_sorted2.csv');
T1 = T([contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1],:);
%%
params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35; % frame sampling rate
params.gsmooth = 0;
params.epochL = 500;
params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
params.padding = 2*params.halfpadding;
rate = 1;
%%
roi_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\roi';
% for kk = 1:size(T1,1)
for kk = 28
    clearvars -except T1 kk params roi_folder rate
    ops = get_session_info(T1,kk);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    serverRoot = expPath(ops.mn, ops.td, ops.en);  
    %% SVD, plot example trace, overlay with pupil and stim time
    [U,V,t,mimg] = get_wf_svd1(serverRoot);
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    nt = numel(t);
    frameRange = 1:params.epochL:nt;
    %%
    U1 = U(1:params.downscale:end,1:params.downscale:end,1:50);
    mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
    mimg2 = padZeros(mimg1,params.halfpadding);
    %%
    fname1 = fullfile(roi_folder,[fname '_roi']);
    load([fname1 '.mat']);
    BW = createMask(roi,mimg2);
    BW2 = BW(121:end-120,121:end-120);
    BW3 = not(BW2);
    %%
    scale = 4;
    U1 = U1(1:scale:end,1:scale:end,:);
    mimg1 = mimg1(1:scale:end,1:scale:end);
    BW3 = BW3(1:scale:end,1:scale:end);
    traceAmp = [];
    for kkk = 1:numel(frameRange)-1
        tic
        if kkk == 1
            frameStart = 1;
            frameTemp = 1:params.epochL;
            dV1 = dV(1:50,frameTemp);
            t1 = t(frameTemp);
            [~,traceAmp1,~] = spiralPhaseMap4(U1,dV1,t1,params,rate,mimg1);
        elseif kkk == numel(frameRange)-1
            frameStart = frameRange(kkk);
            frameTemp = frameRange(kkk):nt;
            dV1 = dV(1:50,frameTemp);
            t1 = t(frameTemp);
            [~,traceAmp1,~] = spiralPhaseMap4(U1,dV1,t1,params,rate,mimg1);
        else
            frameStart = frameRange(kkk);
            frameEnd = frameStart+params.epochL-1;
            if frameEnd+35>nt
                frameTemp = frameStart:frameEnd; % extra 2*35 frames before filter data 
                dV1 = dV(1:50,frameTemp);
                t1 = t(frameTemp);
                [~,traceAmp1,~] = spiralPhaseMap4(U1,dV1,t1,params,rate,mimg1);               
            else
                frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
                dV1 = dV(1:50,frameTemp);
                t1 = t(frameTemp);
                [~,traceAmp1,~] = spiralPhaseMap4(U1,dV1,t1,params,rate,mimg1);
                traceAmp1 = traceAmp1(:,:,1+35:end-35); % reduce 2*35 frames after filter data   
            end
        end
        traceAmp2 = reshape(traceAmp1,size(traceAmp1,1)*size(traceAmp1,2),size(traceAmp1,3));
        traceAmp2(BW3(:),:) = nan;
        traceAmp2_mean = mean(traceAmp2,1,"omitnan");
        traceAmp = [traceAmp,traceAmp2_mean];
        fprintf('Frame %g/%g; time elapsed %g seconds \n', [frameStart,nt, toc])
    end
    save([fname '_amp.mat'],'traceAmp');
end
