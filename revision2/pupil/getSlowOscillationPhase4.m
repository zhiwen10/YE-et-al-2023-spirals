function getSlowOscillationPhase3(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));
BW = logical(projectedAtlas1);                                             % atlas brain boundry binary mask 
[row,col] = find(BW);
brain_index = [col,row];                                                   % position idenx within the brain boundry
%% set params
hist_bin = 40;                                                             % bin size (pixels) for spiral density estimation
pixSize = 0.01;                                                            % pixel resolution mm/pix
pixArea = pixSize^2;                                                       % 2d-pixel resolution mm^2/pix^2
scale  = 1;                                                                % no need to scale for atlas outline here
%%
pixel(1,:) = [590,750]; % SSp-ul
pixel(2,:) = [520,850]; % SSp-ll
pixel(3,:) = [480,960]; % SSp-m
pixel(4,:) = [550,960]; % SSp-n
pixel(5,:) = [682,950]; % SSp-bfd
%% load session table

T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
T = T(T.face &T.eye_DLC,:);
% T = T(logical(T.eye_DLC),:);
% T = T(logical(T.face),:);
% T = T([1,2,6:13],:);
%%
n_bins = 18; % Number of phase bins (20-degree bins)
phase_bins = linspace(-pi, pi, n_bins+1);
spirals_sort =  cell(size(T,1),n_bins);
for kk = 1:size(T,1)                                                 % plot 6 session in one figure, to avoid crowding
    %% session info
    clear spiralsT filteredSpirals trace
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    serverRoot = expPath(mn, td, en);
    [U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    sizeTemplate = [1320,1140];
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    Utransformed1 = Utransformed./mimgtransformed;
    %%
    load(fullfile(data_folder,'spirals','spirals_index',...
    [fname '_motion_energy.mat']));
    image_energy2(isnan(image_energy2)) = 0;
    if numel(image_energy2)<size(t,1)
        image_energy2(numel(image_energy2)+1:size(t,1)) = 0;
    elseif numel(image_energy2)>size(t,1)
        image_energy2 = image_energy2(1:numel(t));
    end 
    trace = image_energy2;
    %%
%     tFile1 = fullfile(serverRoot, 'pupilsize.raw.npy'); 
%     pupil_mean = readNPY(tFile1); 
%     pupil_mean2 = pupil_mean(1:2:end);
%     if numel(pupil_mean2)<size(t,1)
%         pupil_mean2(numel(pupil_mean2)+1:size(t,1)) = 0;
%     elseif numel(pupil_mean2)>size(t,1)
%         pupil_mean2 = pupil_mean2(1:numel(t));
%     end 
%     trace = pupil_mean2;
%%
    freq2 = [0.05,0.5];
    Fs = 35;
    [f1,f2] = butter(4, freq2/(Fs/2), 'bandpass');
    meanTrace_low = filtfilt(f1,f2,trace);
    traceHilbert_low =hilbert(meanTrace_low);
    tracePhase_low = angle(traceHilbert_low);
    traceAmp_low = abs(traceHilbert_low);
    %% only use spiral sequences with at least 2 consecutive frames
    % frame_all = numel(t);    
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);                                          % sprial length > = 2 frames
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));                  % transform sprials to atlas space
    filteredSpirals(:,1:2) = round(spiralsT);    
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % find spirals within the brain boundry
    filteredSpirals = filteredSpirals(lia,:);   
    %%control_all1 = squeeze(count_sample(:,pairs2(i,:)));  
    for mm = 1:n_bins
        clear index a b       
        index = find(tracePhase_low >= phase_bins(mm) & tracePhase_low < phase_bins(mm+1));
        [a,b] = ismember(filteredSpirals(:,5),index);
        spirals_sort{kk,mm} = filteredSpirals(a,:);
    end  
end
%%
save('spirals_sort_by_slow_phase3.mat','spirals_sort');