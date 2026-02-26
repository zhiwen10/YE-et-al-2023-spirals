function getSpiralsWhiskerMeanMaps(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
BW2 = BW(1:2:end,1:2:end);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
%% right SSp index
clear areaPath
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
hemi = 'right';
Utransformed = projectedAtlas1;
scale = 1;
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp2 = [col,row];
%%
spiral_folder = fullfile(data_folder,'whisker','spirals_all','sprials_grouping');
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
load(fullfile(data_folder,'whisker','whisker_mean_maps','whisker_spirals_mean_all.mat'));
wf_mean3 = imresize(wf_mean2,[660,570]);
%% tracePhase
wf1 = reshape(wf_mean3,size(wf_mean3,1)*size(wf_mean3,2),size(wf_mean3,3));
meanTrace = wf1 -mean(wf1 ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
meanTrace2 = meanTrace';
meanTrace2 = reshape(meanTrace2,size(wf_mean3,1),size(wf_mean3,2),size(wf_mean3,3));
traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = tracePhase';
tracePhase = reshape(tracePhase, size(wf_mean3,1),size(wf_mean3,2),size(wf_mean3,3));
%% spiral detection
freq = [2,8];                                                              % data filtering frequency range
rate = 1;                                                                  % set to 1, if no upsampling in time
params.downscale = 1;                                                      % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0;                                                        % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35;                                                            % frame sampling rate
params.halfpadding = 120;                                                  % half padding should be the max spiral detection radius 
params.padding = 2*params.halfpadding;
params.th = 1:36:360;                                                      % 10 evenly-spaced angles in a circle
params.rs = 10:5:20;                                                       % candidate spiral check at small radii [10, 15, 20 pixels, 17.3um/pixel]
params.gridsize = 10;                                                      % search grid resolution at 10 pixels/grid
params.spiralRange = linspace(-pi,pi,5);                                   % evenly-spaced 4 phase quadrants of a circle
params.gsmooth = 0;                                                        % no spatial smoothing 
params.epochL = 141;                                                       % frames per processing batch 
params.nt = size(tracePhase,3);                                            % total frames
params.frameRange = 1:params.nt;                                           % frameRange excludeing 1 second edges
params.dThreshold = 15;                                                    % distance threshold for sprial center clusting
params.rsRCheck = 10:10:100;                                               % final spiral radius check range
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);        % roi grids for refined search
%% generate coarse search grids with zeros padded at the edges
params.xsize = size(tracePhase,1);
params.ysize = size(tracePhase,2);
xsizePadded = params.xsize+params.padding; 
ysizePadded = params.ysize+params.padding;
[xx,yy] = meshgrid(...
    min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
    min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
params.xx = xx; params.yy = yy;
%% draw brain mask roi
mimg1 = squeeze(tracePhase(:,:,100));
mimg2 = padZeros(mimg1,params.halfpadding);
figure; 
ax1 = imagesc(mimg2);
colormap(colorcet('C06'));
roi = drawpolygon;
%%
tf = inROI(roi,params.xx(:),params.yy(:));
params.xxRoi = params.xx(tf);                                                % only use the grids that inside the roi to save time
params.yyRoi = params.yy(tf);                                                % only use the grids that inside the roi to save time         
%% spiral detection
pwAll = [];
pwAll1 = []; pwAll2 = []; 
pwAll3 = []; pwAll4 = []; pwAll5 = [];
tracePhase1 = padZeros(tracePhase,params.halfpadding);                       % pad tracephase with edge zeros
nframe = size(tracePhase1,3);
frameStart = 1;
for frame = 1:nframe
    A = squeeze(tracePhase1(:,:,frame));
    pwAll1 = spiralAlgorithm(A,params);                                      % coarse search
    frameN = frame+frameStart-1;
    pwAll2 = checkClusterXY(pwAll1,params.dThreshold);                      % spatial clustering of nearby duplicate spiral centers
    pwAll3 = doubleCheckSpiralsAlgorithm(A,pwAll2,params);                  % double check the mean cluster centers from the last step are still spiral centers   
    pwAll4 = spatialRefine(A,pwAll3,params);                                % refined spiral search based on candidate sprial centers
    [pwAll5] = spiralRadiusCheck2(A,pwAll4,params);                         % find spiral radius and traveling direction
    if not(isempty(pwAll5))
        pwAll5(:,end+1) = frameN;                                            % attach frame ID label
    end
    pwAll = [pwAll;pwAll5];                                                  % concatenate all sprials from the current batch
end
pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;                            % recalculate spiral 2d coordinates without padding   
filteredSpirals = pwAll(pwAll(:,3)>=40,:);                                   % only use sprials with radius >40 pixels, based on 3d-fft
filteredSpirals =unique(filteredSpirals, 'rows');                            % get rid of duplication, in case any
filteredSpirals = sortrows(filteredSpirals,5);                               % sort based on frame number
[archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);            % main grouping algorithm
%%
save(fullfile(save_folder,'whisker_spirals_group_fftn.mat'),'archiveCell');
save(fullfile(save_folder,'whisker_evoked_map_roi.mat'),'roi'); 