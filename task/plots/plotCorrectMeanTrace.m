function h5c = plotCorrectMeanTrace(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%% load data from an example an session
mn = 'ZYE_0091';
T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
T1 = T_session(T_session.label == "task",:);
T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
sessions = size(T1,1);
kk = 1;
mn = char(T1.MouseID{kk});
tda = T1.date(kk);  
en = T1.folder(kk); 
block_en = T1.block_folder(kk); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
%% load atlas registration
downscale = 8;
load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(V(1:50,:));
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 0.1;
pixel(1,:) = [870,850]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [660,950]; % SSp-bfd
pixel = round(pixel/params.downscale);
%%
Utransformed = Utransformed(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
%%
load(fullfile(data_folder,'task','task_mean_maps','task_mean_maps_all_mice.mat'));
%%
t1 = -1:1/35:1;
tq = 1:rate:numel(t1);
qt = interp1(1:numel(t1),t1,tq);
%%
color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
frames_to_plot = 18;
first_frame =36;
last_frame = first_frame+frames_to_plot;
t1a = t1(first_frame); t1b = t1(last_frame);
first_frame1 = find(qt-t1a>0, 1, 'first');
last_frame1 = find(qt-t1b>0, 1, 'first');

example_frame_to_plot = first_frame+2;
t1c = t1(example_frame_to_plot); 
frame = find(qt-t1c>0, 1, 'first');
%%
scale3 = 5/8;
load(fullfile(data_folder,'tables','task_mask_all_mice.mat'));
BW2 = BW2(1:downscale:end,1:downscale:end);
%%
h5c = figure('Renderer', 'painters', 'Position', [100 100 400 300]);
ax1 = subplot(1,2,1);
low = prctile(mimgtransformedRGB(:),2);
high = prctile(mimgtransformedRGB(:),98);
im1= imshow(mimgtransformedRGB,[low high]);
set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
% cb1 = colorbar(ax1);
% pos1 = cb1.Position;
% cb1.Position = [pos1(1)+0.1, pos1(2),pos1(3),pos1(4)/2];
hold on;
scatter(ax1,pixel(:,2),pixel(:,1),16,color2(3:end,:),'filled');
%
for kk = 1
    % -2:2 seconds 141 samples
    trace_mean_current = squeeze(trace_mean_all(:,:,:,kk));
    trace_mean_current = imresize(trace_mean_current,[165, 143]);
    freq = [2,8];
    [trace_mean3,traceFilt3,tracePhase3] = taskTrace_upsample(trace_mean_current,freq,rate);
    % only use -1:1 seconds 70 samples
    trace2d = trace_mean3(:,:,1+35/rate:end-35/rate); 
    traceFilt2d = traceFilt3(:,:,1+35/rate:end-35/rate); 
    tracePhase1 = tracePhase3(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
    tracePhase1(isnan(tracePhase1(:))) = 0;
    %
    trace_filt = zeros(3,size(trace2d,3));
    trace_raw = zeros(3,size(trace2d,3));
    trace_phase = zeros(3,size(trace2d,3));
    for i = 1:7
        trace_raw(i,:) = trace2d(pixel(i,1),pixel(i,2),:);
        trace_filt(i,:) = traceFilt2d(pixel(i,1),pixel(i,2),:);
        trace_phase(i,:) = tracePhase1(pixel(i,1),pixel(i,2),:);
    end
    %%
    ax2 = subplot(1,2,kk+1);
    scalea = 1.5;
    for i = 1:7
        % plot(ax2,qt(280:561),trace_raw(i,280:561)*scalea*100+2*(i-1),'k','lineWidth',1);
        plot(ax2,qt,trace_raw(i,:)*scalea*100+2*(i-1),'k','lineWidth',1);
        hold on; 
        [pks,locs] = findpeaks(trace_raw(i,:),'MinPeakProminence',0.001);
        scatter(qt(locs),trace_raw(i,locs)*scalea*100+2*(i-1),6,'r','filled');
    end
    ax2.FontSize = 8; 
    xlim([qt(280),qt(561)]);

    set(ax2,'YTickLabel',[]);
    for i = 1:7
        text(qt(311)-0.1,(i-1)*2,nameList{i},'Color', color2(i+2,:),'Interpreter','None');
    end
    xlabel('Time (s)','fontSize',9);
    line_ax1 = xline(qt(first_frame1),'--');
    line_ax2 = xline(qt(last_frame1),'--');
    hold on;
    plot([qt(280)+0.5,qt(280)+0.5],[0,scalea],'r');
    xlim([-0.2,0.6]);
    xticks([-0.2,0,0.2,0.4,0.6]);
    xticklabels({'-0.2','0','0.2','0.4','0.6'});
    hold on;
    plot([0,0.2],[0,0],'k','lineWidth',2);
end
%%
print(h5c,fullfile(save_folder,'Fig5c_task_correct_mean_trace2.pdf'),'-dpdf', '-bestfit', '-painters');