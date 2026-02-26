githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
%% load widefield data
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % read the session table
kk = 1;                                                                    % take a look at first session in the table
mn = T.MouseID{kk};                                                        % mouse name
tda = T.date(kk);                                                          % session date
en = T.folder(kk);                                                         % experiment folder
tdb = datestr(tda,'yyyymmdd');                                             % session date string
subfolder = [mn '_' tdb '_' num2str(en)];                                  % full name identifier for the session
session_root = fullfile(data_folder,'spirals','svd',subfolder);              % svd data location for the session
[U,V,t,mimg] = loadUVt1(session_root);                                     % load svd components (U,V,t) from svd data folder
dV = [zeros(size(V,1),1) diff(V,[],2)];                                    % take derivative of the V components
U1 = reshape(U(:,:,1:50),size(U,1)*size(U,2),50);                          % only use 50 components (>99% variance), and reshape to 2d matrix
frames = [200:208];                                                        % let's check random 8 frames
trace = U1*dV(1:50,frames);                                                % reconstruct svd data to image space 
trace = reshape(trace,size(U,1),size(U,2),size(trace,2));                  % reshape images to x*y*t matrix
trace = trace./mimg;                                                       % normalize by mean intensity, to get df/f
frameN = numel(frames);
cmin = min(trace(:)); cmax = max(trace(:));
fig = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);      
for i = 1:frameN
    ax(i) = subplot(1,frameN,i);
    imagesc(squeeze(trace(:,:,i)));                                        % plot example random 8 frames of raw data
    caxis([cmin,cmax]);                                                    % same color scale for all images
    axis image; axis off;   
end
h = axes(fig,'visible','off');                                  
c = colorbar(h,'Position',[0.93 0.35 0.01 0.3]);                           % plot color scale bar
caxis([cmin,cmax]);
c.Label.String = 'dF/F';
%% load axon morphology data
load(fullfile(data_folder,'axons','all_cell_with_parents.mat'));             % load singel cell morphology dataset
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % load 10um resolution atlas horizontal view outlines 
id = [82,806,799];                                                         % let's take a look at these 3 random cells
somaType = 1;                                                              % soma coordinates are identified by type1 on the first data column
axonType = 2;                                                              % axons coordinates are identified by type2 on the first data column
scale = 1;
fig = figure('Renderer', 'painters', 'Position', [100 100 1000 300]);   
for i = 1:3
    ce1l1 = allCoords{id(i)};
    axonCoords = ce1l1(ce1l1(:,1)== axonType,2:4);                         % find 3d coordinates of all axon points
    somaCoords = ce1l1(ce1l1(:,1)== somaType,2:4);                         % find 3d coordinates of all soma points
    subplot(1,3,i);
    overlayOutlines(coords,scale,[0.8,0.8,0.8]);                           % plot brain atlas outlines in horizontal view
    set(gca, 'YDir','reverse');
    plot(axonCoords(:,3)/10, axonCoords(:,1)/10,...
        '.', 'Color','k','MarkerSize',1);                                  % plot all axons in 2d horizontal view, for simplicity
    hold on;
    plot(somaCoords(:,3)/10, somaCoords(:,1)/10,...                        % plot soma location
        '.', 'Color','r','MarkerSize',12); 
    axis image; axis off;
end

%% load simultaneous widefield and ephys data 
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions.csv'));
kk = 1;
ops = get_session_info2(T,kk,data_folder);                                 % get session info in a structure array
[U,V,t,mimg] = loadUVt1(ops.session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                    % get derivative of V
[sp] = loadKSdir2(ops.session_root);                                       % load spike sorting results
[syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t);                     % synchronization step, change widefield t to ephys spiketime t
WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));                                % discard nan
nVc = 50;                                                                  % how many components to use, 50 components contain >99% variance
dV1 = double(dV(1:nVc,~isnan(WF2ephysT1)));  
V1 = double(V(1:nVc,~isnan(WF2ephysT1)));
tt1 = WF2ephysT1(~isnan(WF2ephysT1));
[MUA_std] = get_MUA_bin(sp,WF2ephysT);                                     % bin spiking data based on widefield T and  sampling rate
dV1 = double(dV1(1:nVc,~isnan(WF2ephysT1)));                               % now dV1 and MUA_std should have the same length in time dimension
fig = figure('Renderer', 'painters', 'Position', [100 100 1000 400]); 
pixelA = [400,250];                                                        % let's look at an example pixel in the left RSP, with 2-8Hz oscillation      
traceA = squeeze(U(pixelA(1),pixelA(2),1:nVc))'*dV1;                        % reconstruct pixel trace from svd 
traceA = traceA./mimg(pixelA(1),pixelA(2));                                % normalize by mean intensity, to get df/f
MUA_sum = sum(MUA_std,1);                                                  % sum all MUA activity
MUA_sum = MUA_sum./max(MUA_sum);                                           % normalize by max activity
subplot(1,3,1);
imagesc(mimg);
axis image; axis off;
hold on;
scatter(pixelA(2),pixelA(1),12,'r','filled');                              % plot the RSP pixel on the widefield image
subplot(1,3,[2,3]);
plot(tt1,traceA*100,'r');                                                  % plot widefield trace
hold on;
plot(tt1,MUA_sum+1.5,'k');                                                 % plot spiking MUA trace
xlim([390,410]);                                                           % let's zoom in on 20 seconds
legend({'widefield','MUA'})
xlabel('Time (s)');