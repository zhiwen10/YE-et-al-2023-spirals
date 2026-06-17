function getSpiralDetectionCurl(T, data_folder, save_folder)
% GETSPIRALDETECTIONCURL  Detect cortical spiral waves using multi-scale curl.
%
% ALGORITHM OVERVIEW:
%   Replaces the original grid-search / winding-number approach with a
%   fully vectorised curl-of-phase-gradient method.
%
%   For each frame's bandpass-filtered Hilbert phase map the pipeline is:
%
%   (1) MULTI-SCALE CURL MAP
%       For each stencil radius r in rs_curl, compute the discrete curl at
%       every pixel simultaneously using the 4-corner formula derived from
%       Green's theorem:
%
%         curl_r(y,x) =  wrap[ phi(y+r, x+r) - phi(y-r, x+r) ]   % dPhi_dy at x+r
%                      - wrap[ phi(y+r, x-r) - phi(y-r, x-r) ]   % dPhi_dy at x-r
%                      - wrap[ phi(y+r, x+r) - phi(y+r, x-r) ]   % dPhi_dx at y+r
%                      + wrap[ phi(y-r, x+r) - phi(y-r, x-r) ]   % dPhi_dx at y-r
%
%       All wrapping is done via angle(exp(1i * ...)) for numerical
%       stability near the +/-pi discontinuity.  Each term is a finite
%       difference of the wrapped phase gradient; together they approximate
%       d(dPhi/dy)/dx - d(dPhi/dx)/dy, the z-component of the curl.
%
%       By Stokes' theorem the line integral of grad(phi) around any closed
%       loop encircling a phase singularity equals +/-2*pi (for winding
%       number +/-1).  The 4-corner formula computes exactly this
%       circulation per unit area.  Verified: for phi = atan2(y,x) (CCW)
%       the formula gives +2*pi at the origin and 0 elsewhere.
%
%   (2) MULTI-SCALE VOTE MAP (noise robustness)
%       A true spiral centre produces |curl_r| ~ 2*pi at EVERY scale r
%       because the 2*pi phase winding is a topological invariant.
%       Noise pixels produce large curl only at small r (where isolated bad
%       pixels dominate the 4-corner differences) but not at large r where
%       local fluctuations average out.
%
%       Vote rule per scale r:
%         votePos += 1  if  curl_r > +curl_thresh   (CCW evidence)
%         voteNeg += 1  if  curl_r < -curl_thresh   (CW  evidence)
%
%       A pixel is a candidate centre only if max(votePos, voteNeg) >=
%       min_votes AND the votes are all in the same direction (sign
%       consistency).  This "K-out-of-N scales" criterion mirrors the
%       original "2 out of 3 radii" criterion.
%
%   (3) CANDIDATE DETECTION
%       Apply brain ROI mask (eroded by max stencil radius to avoid border
%       artefacts from zero-padding).  Find local maxima of the vote map
%       using morphological non-maximum suppression.
%
%   (4) CLUSTERING + RADIUS / DIRECTION (reuse original helpers)
%       Cluster nearby maxima with checkClusterXY.  Determine spiral radius
%       and rotation direction with spiralRadiusCheck2, with the direction
%       confirmed / overridden from the curl sign map.
%
% INPUTS:
%   T            - session table with columns: MouseID (cell of strings),
%                  date (datetime or datenum), folder (numeric)
%   data_folder  - root path; must contain:
%                    <data_folder>/spirals/svd/<session>/   (SVD files)
%                    <data_folder>/spirals/full_roi/         (ROI .mat files)
%   save_folder  - destination folder for output CSV files
%
% OUTPUT CSV columns (one row per detected spiral):
%   spiral_center_x    - x pixel coordinate (un-padded image space)
%   spiral_center_y    - y pixel coordinate (un-padded image space)
%   spiral_radius      - estimated outer radius (pixels)
%   spiral_direction   - rotation: +1 = CCW, -1 = CW
%   spiral_frame       - frame index in the original recording

for kk = 1:size(T, 1)

    %% session identifiers
    mn  = T.MouseID{kk};
    tda = T.date(kk);
    en  = T.folder(kk);
    tdb = datestr(tda, 'yyyymmdd');

    %% load SVD components
    subfolder    = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder, 'spirals', 'svd', subfolder);
    [U, V, t, ~] = loadUVt1(session_root);
    dV = [zeros(size(V,1), 1)  diff(V, [], 2)];

    %% set detection parameters
    freq   = [2, 8];   % bandpass filter range (Hz)
    rate   = 1;        % no temporal upsampling
    params = setSpiralDetectionParamsCurl(U, t);

    %% build binary ROI mask over the padded image domain
    fname1  = [mn '_' tdb '_' num2str(en) '_roi'];
    roiFile = load(fullfile(data_folder, 'spirals', 'full_roi', [fname1 '.mat']));
    roi     = roiFile.roi;
    U1      = U(1:params.downscale:end, 1:params.downscale:end, 1:50);
    params.roiMask = buildRoiMask(roi, params, size(U1,1), size(U1,2));

    %% run curl-based spiral detection
    pwAll = spiralDetectionCurlAlgorithm(U1, dV, t, params, freq, rate);

    %% write results as CSV
    fname = [mn '_' tdb '_' num2str(en) '_spirals_curl.csv'];
    if isempty(pwAll)
        pwAll = zeros(0, 5);
    end
    T1 = array2table(pwAll, 'VariableNames', ...
        {'spiral_center_x', 'spiral_center_y', ...
         'spiral_radius', 'spiral_direction', 'spiral_frame'});
    writetable(T1, fullfile(save_folder, fname));
    fprintf('Session %s: %d spirals detected\n', subfolder, size(pwAll, 1));

end
end


% =========================================================================
%  PARAMETER SETUP
% =========================================================================

function params = setSpiralDetectionParamsCurl(U, t)
% Detection parameters.  Inherits timing / padding / filter settings from
% the original pipeline; adds curl-specific parameters.

% --- shared with original pipeline ---
params.downscale   = 1;                     % spatial downscale factor
params.lowpass     = 0;                     % 0 = bandpass, 1 = lowpass
params.Fs          = 35;                    % camera frame rate (Hz)
params.halfpadding = 120;                   % zero-pad half-width (pixels)
params.padding     = 2 * params.halfpadding;
params.gsmooth     = 0;                     % no spatial smoothing
params.epochL      = 100;                   % frames per processing batch
params.nt          = numel(t);
params.frameN1     = round((params.nt - 70) / 100) * 100;
params.frameRange  = 36 : params.epochL : (36 + params.frameN1 - 1);

% --- kept for spiralRadiusCheck2 (radius & direction) ---
params.th          = 1:36:360;              % sample angles on circle
params.spiralRange = linspace(-pi, pi, 5); % phase quadrant edges
params.rsRCheck    = 10:10:100;             % radius sweep range (pixels)
params.dThreshold  = 15;                    % spatial clustering distance (pixels)

% --- curl-specific parameters ---
% Stencil radii: small r catches fine structure but is noise-sensitive;
% large r is robust to noise but may miss spirals smaller than 2r.
% Range chosen to match original rs=[10,15,20] search radii.
params.rs_curl     = [5, 10, 15, 20, 25];  % stencil half-widths (pixels)

% Minimum |curl| per scale to count as a vote.
% For a true spiral: curl ~ 2*pi ~ 6.28 at every scale.
% For noise: curl ~ 0 on average.  Threshold at pi is conservative.
params.curl_thresh = pi;                    % radians (raw 4-corner formula units)

% Require at least this many scales to agree in sign.
% 3 out of 5 mirrors the original "2 out of 3 radii" criterion.
params.min_votes   = 3;                     % minimum consistent-sign votes

% Non-maximum suppression disc radius.
% Should be ~= original gridsize (10 px) to yield one candidate per spiral.
params.nms_radius  = 10;                    % pixels

U1 = U(1:params.downscale:end, 1:params.downscale:end, 1:50);
params.xsize = size(U1, 1);
params.ysize = size(U1, 2);
end


% =========================================================================
%  ROI MASK
% =========================================================================

function roiMask = buildRoiMask(roi, params, xsize, ysize)
% Build a logical mask over the full padded image.
% Pixels are true when they are inside the brain ROI AND at least
% max(rs_curl) pixels from the padded edge (to prevent border artefacts
% caused by circshift wrapping into the zero-padded region).

xsizePad = xsize + params.padding;   % rows (y direction)
ysizePad = ysize + params.padding;   % cols (x direction)

% Build coordinate grids in padded image space.
% meshgrid convention: xx(r,c) = c (column = x), yy(r,c) = r (row = y).
[xx_pad, yy_pad] = meshgrid(1:ysizePad, 1:xsizePad);

% inROI expects (roi, x, y) = (roi, col, row) — same as original pipeline.
tf      = inROI(roi, xx_pad(:), yy_pad(:));
roiMask = reshape(tf, xsizePad, ysizePad);

% Erode by the largest stencil radius so the 4-corner formula never reads
% from the zero-padded region for any candidate centre inside the mask.
se      = strel('disk', max(params.rs_curl));
roiMask = imerode(roiMask, se);
end


% =========================================================================
%  MAIN DETECTION LOOP
% =========================================================================

function pwAll = spiralDetectionCurlAlgorithm(U1, dV, t, params, freq, rate)
% Iterate over frame batches, compute phase maps, and run per-frame curl
% detection.  Mirrors the structure of spiralDetectionAlgorithm.m.

pwAll = [];

for kkk = 1:numel(params.frameRange) - 1
    tic

    frameStart = params.frameRange(kkk);
    frameEnd   = frameStart + params.epochL - 1;
    frameTemp  = (frameStart - 35) : (frameEnd + 35);   % extra frames for filter edges

    dV1 = dV(1:50, frameTemp);
    t1  = t(frameTemp);

    % Bandpass filter + Hilbert transform → phase map (reuse original function)
    [~, ~, tracePhase1] = spiralPhaseMap_freq(U1, dV1, t1, params, freq, rate);
    tracePhase1 = tracePhase1(:,:, 1+35:end-35);        % trim filter edges
    tracePhase  = padZeros(tracePhase1, params.halfpadding);
    nframe      = size(tracePhase, 3);

    for frame = 1:nframe
        frameN = frame + frameStart - 1;
        phi    = squeeze(tracePhase(:,:, frame));        % 2-D phase map, radians in [-pi, pi]

        % --- Step 1: compute multi-scale curl vote map ---
        [voteMap, dirMap] = computeCurlVoteMap(phi, params);

        % --- Step 2: apply brain ROI mask ---
        voteMap(~params.roiMask) = 0;

        % --- Step 3: find local maxima with sufficient votes ---
        candidates = findCurlCandidates(voteMap, params);
        if isempty(candidates), continue, end

        % --- Step 4: cluster nearby candidates (reuse original helper) ---
        candidates = checkClusterXY(candidates, params.dThreshold);
        if isempty(candidates), continue, end

        % --- Step 5: determine radius and direction (reuse original helper) ---
        spirals = spiralRadiusCheck2(phi, candidates, params);
        if isempty(spirals), continue, end

        % Override/confirm direction from curl sign map (more direct than
        % the winding-integral sign used in spiralRadiusCheck2).
        for is = 1:size(spirals, 1)
            px = spirals(is, 1);   % column (x)
            py = spirals(is, 2);   % row    (y)
            if py >= 1 && py <= size(dirMap,1) && ...
               px >= 1 && px <= size(dirMap,2)
                d = dirMap(py, px);
                if d ~= 0
                    spirals(is, 4) = d;   % +1 CCW, -1 CW
                end
            end
        end

        spirals(:, end+1) = frameN;                      %#ok<AGROW>
        pwAll = [pwAll; spirals];                        %#ok<AGROW>
    end

    fprintf('Frame %g / %g;  elapsed %.1f s\n', frameStart, params.frameN1, toc);
end

% Convert padded pixel coordinates back to original image coordinates
if ~isempty(pwAll)
    pwAll(:,1) = pwAll(:,1) - params.halfpadding;
    pwAll(:,2) = pwAll(:,2) - params.halfpadding;
end
end


% =========================================================================
%  CORE ALGORITHM: MULTI-SCALE CURL VOTE MAP
% =========================================================================

function [voteMap, dirMap] = computeCurlVoteMap(phi, params)
% COMPUTECURLVOTEMAP  Vectorised multi-scale curl computation.
%
% For each stencil radius r, samples the phase at the 4 corners of a
% (2r x 2r) square centred at every pixel and computes the discrete curl:
%
%   Let:  A = phi(y+r, x+r)   [circshift(phi, [-r, -r])]
%         B = phi(y+r, x-r)   [circshift(phi, [-r, +r])]
%         C = phi(y-r, x+r)   [circshift(phi, [+r, -r])]
%         D = phi(y-r, x-r)   [circshift(phi, [+r, +r])]
%
%   curl_r = [wrap(A-C) - wrap(B-D)] - [wrap(A-B) - wrap(C-D)]
%          = [dPhi_dy(x+r) - dPhi_dy(x-r)] - [dPhi_dx(y+r) - dPhi_dx(y-r)]
%
%   where wrap(*) = angle(exp(1i * *))  keeps values in (-pi, pi].
%
%   Verification (analytical):
%     phi = atan2(y,x) (CCW spiral at origin):
%       A=pi/4, B=3pi/4, C=-pi/4, D=-3pi/4  (at stencil r, centre pixel)
%       curl = [pi/2 - (-pi/2)] - [-pi/2 - pi/2] = pi + pi = +2*pi ✓
%     phi = -atan2(y,x) (CW spiral):
%       curl = -2*pi ✓
%     Plane wave phi = k*x: all differences cancel, curl = 0 ✓
%
%   NOISE ROBUSTNESS:
%     A noisy pixel perturbs only 1-2 corners of the 4-corner stencil.
%     At small r the perturbation is large relative to the curl signal.
%     At large r the 4 corners are far from the noise source and the curl
%     returns to 0.  A true spiral centre gives |curl| ~ 2*pi at ALL r.
%     The vote map captures this: only pixels with K/N consistent votes
%     (same sign, above threshold) are kept as candidates.
%
% OUTPUTS:
%   voteMap  - uint8 map, range [0, numel(rs_curl)]
%              value = number of scales with consistent-sign curl above thresh
%   dirMap   - int8 map: +1 CCW, -1 CW, 0 ambiguous

rs_curl   = params.rs_curl;
threshold = params.curl_thresh;

votePos = zeros(size(phi), 'uint8');   % CCW votes per pixel
voteNeg = zeros(size(phi), 'uint8');   % CW  votes per pixel

for r = rs_curl
    % 4 corner phase maps via circular shift.
    % MATLAB circshift convention: circshift(A, [dr, dc]) shifts rows DOWN
    % by dr and cols RIGHT by dc.  A negative shift looks ahead:
    %   circshift(phi, [-r,  0]) at (i,j) gives phi(i+r, j)   → phi(y+r, x)
    %   circshift(phi, [+r,  0]) at (i,j) gives phi(i-r, j)   → phi(y-r, x)
    %   circshift(phi, [ 0, -r]) at (i,j) gives phi(i, j+r)   → phi(y, x+r)
    %   circshift(phi, [ 0, +r]) at (i,j) gives phi(i, j-r)   → phi(y, x-r)
    A = circshift(phi, [-r, -r]);   % phi(y+r, x+r)
    B = circshift(phi, [-r, +r]);   % phi(y+r, x-r)
    C = circshift(phi, [+r, -r]);   % phi(y-r, x+r)
    D = circshift(phi, [+r, +r]);   % phi(y-r, x-r)

    % Wrapped differences via complex exponential (avoids explicit modular
    % arithmetic; numerically stable at the +/-pi discontinuity).
    w = @(p, q) angle(exp(1i .* (p - q)));

    % Discrete curl = d(dPhi_dy)/dx - d(dPhi_dx)/dy
    %   d(dPhi_dy)/dx  ∝  wrap(A-C) - wrap(B-D)   [gradient of dPhi/dy in x]
    %   d(dPhi_dx)/dy  ∝  wrap(A-B) - wrap(C-D)   [gradient of dPhi/dx in y]
    curl_r = (w(A,C) - w(B,D)) - (w(A,B) - w(C,D));

    votePos = votePos + uint8(curl_r >  threshold);
    voteNeg = voteNeg + uint8(curl_r < -threshold);
end

% voteMap: how many scales agree (whichever direction has more votes wins)
voteMap = max(votePos, voteNeg);

% dirMap: +1 CCW, -1 CW, 0 if tied or zero
dirMap = int8(sign(int16(votePos) - int16(voteNeg)));
end


% =========================================================================
%  CANDIDATE LOCALISATION: NON-MAXIMUM SUPPRESSION
% =========================================================================

function candidates = findCurlCandidates(voteMap, params)
% Find local maxima of voteMap that meet the minimum vote threshold.
%
% Uses morphological dilation for non-maximum suppression (NMS):
% a pixel is a local maximum if its vote count equals the maximum within
% a disc of radius nms_radius.  This collapses the broad curl blob
% produced by a single spiral (width ~ sqrt(2)*r) down to its peak,
% giving a single candidate per spiral.
%
% Returns candidates as [x, y] = [col, row] pairs, matching the
% convention used by checkClusterXY and spiralRadiusCheck2.

candidates = [];
minVotes   = params.min_votes;
nmsR       = params.nms_radius;

if ~any(voteMap(:) >= minVotes)
    return
end

% Morphological NMS: pixel is max iff it equals the dilated image.
se       = strel('disk', nmsR);
localMax = (imdilate(double(voteMap), se) == double(voteMap)) & ...
           (voteMap >= minVotes);

[rows, cols] = find(localMax);
if isempty(rows)
    return
end

% Return as [x, y] = [col, row]
candidates = [cols, rows];
end
