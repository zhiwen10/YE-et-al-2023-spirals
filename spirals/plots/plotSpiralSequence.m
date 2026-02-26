function h1b = plotSpiralSequence(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%% load data from an example an session
mn = 'ZYE_0012';
td = '2020-10-16';
tdb = datestr(td,'yyyymmdd');
en = 5;
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
load(fullfile(data_folder,'tables','mask_ZYE12.mat'));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'tables',[fname '_tform.mat']));               % load atlas transformation matrix tform;
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
pixel(1,:) = [900,800]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [682,950]; % SSp-bfd
pixel = round(pixel/params.downscale);
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%% hilbert transform
freq = [2,8];
tStart = 1681; tEnd = 1686; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(:,frameTemp);
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap_freq2(Utransformed,dV1,t,params,freq,rate);
trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimgtransformed; 
tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
%%
trace = zeros(7,size(trace2d1,3));
for i = 1:7
    trace(i,:) = trace2d1(pixel(i,1),pixel(i,2),:);
end
%%
t1 = t(frameStart:frameEnd);
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
color1 = cbrewer2('seq','YlOrRd',7);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
%%
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
frames = 18;
col = 9;
first_frame =23;

dff = trace2d1*100;
raw_min = min(dff(:)); raw_max = max(dff(:));
phase_min = -pi; phase_max = pi;

h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
frame_count = 1;
for k = 1:frames
    frame  =first_frame+(k-1);
    ax3 = subplottight(4,col,k);
    im_raw = imagesc(dff(:,:,frame));
    colormap(ax3, parula)
    axis image; axis off;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([raw_min*0.6,raw_max*0.6]);
    text_ax1 = title(['frame',num2str(frame_count)],'fontSize',8);
    
    ax4 = subplottight(4,col,frames+k);
    im_phase = imagesc(tracePhase1(:,:,frame));
    colormap(ax4,colorcet('C06'));
    
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    caxis([phase_min,phase_max]);
    text_ax2 = title(['frame',num2str(frame_count)],'fontSize',8);
    frame_count = frame_count+1;
end
%%
print(h1b, fullfile(save_folder,'Fig1b_example_spiral_sequence.pdf'),...
    '-dpdf', '-bestfit', '-painters');
