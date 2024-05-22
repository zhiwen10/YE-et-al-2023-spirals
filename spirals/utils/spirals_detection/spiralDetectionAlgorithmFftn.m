function [pwAll] = spiralDetectionAlgorithmFftn(U1,dV,t,params,freq,rate) 
pwAll = [];
for kkk = 1:numel(params.frameRange)-1                                     % loop through batchs of frames
    tic
    pwAll1 = []; pwAll2 = []; 
    pwAll3 = []; pwAll4 = []; pwAll5 = [];
    frameStart = params.frameRange(kkk);
    frameEnd = frameStart+params.epochL-1;
    frameTemp = frameStart-35:frameEnd+35;                                 % extra 2*35 frames before filter data 
    dV1 = dV(1:50,frameTemp);
    t1 = t(frameTemp);
    [trace2d1,traceAmp1,tracePhase1] = ...
        spiralPhaseMap_fftn_freq(U1,dV1,t1,params,freq,rate);              % filter data at a pre-defined frequency range
    tracePhase1 = tracePhase1(:,:,1+35:end-35);                            % reduce 2*35 frames after filter data 
    tracePhase = padZeros(tracePhase1,params.halfpadding);                 % pad tracephase with edge zeros
    nframe = size(tracePhase,3);
    for frame = 1:nframe
        A = squeeze(tracePhase(:,:,frame));
        pwAll1 = spiralAlgorithm(A,params);                                % coarse search
        frameN = frame+frameStart-1;
        pwAll2 = checkClusterXY(pwAll1,params.dThreshold);                 % spatial clustering of nearby duplicate spiral centers
        pwAll3 = doubleCheckSpiralsAlgorithm(A,pwAll2,params);             % double check the mean cluster centers from the last step are still spiral centers   
        pwAll4 = spatialRefine(A,pwAll3,params);                           % refined spiral search based on candidate sprial centers
        [pwAll5] = spiralRadiusCheck2(A,pwAll4,params);                    % find spiral radius and traveling direction
        if not(isempty(pwAll5))
            pwAll5(:,end+1) = frameN;                                      % attach frame ID label
        end
        pwAll = [pwAll;pwAll5];                                            % concatenate all sprials from the current batch
    end
    fprintf('Frame %g/%g; time elapsed %g seconds \n',...
        [frameStart,params.frameN1, toc])
end
pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;                            % recalculate spiral 2d coordinates without padding   