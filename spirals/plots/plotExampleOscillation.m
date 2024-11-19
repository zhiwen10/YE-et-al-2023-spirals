function hs1ab = plotExampleOscillation(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
mn = 'ZYE_0052';
td = '2021-12-18';
tdb = datestr(td,'yyyymmdd');
en = 2;
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];
%%
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','spirals_example',[fname '_lick_wheel.mat']));
%% registration
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform;
%%
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
%%
params.downscale = 1;
Utransformed = Ut(1:params.downscale:end,1:params.downscale:end,:);
mimgtransformed = mimgt(1:params.downscale:end,1:params.downscale:end);
mimgtransformedRBG = mimgtransformed/max(mimgtransformed(:));
mimgtransformedRGB = cat(3, mimgtransformedRBG,mimgtransformedRBG,mimgtransformedRBG);
BW = logical(projectedAtlas1);
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
%%
Ut1 = Ut./mimgt;
BW = logical(projectedAtlas1);
%% ZYE52 points
points(1,:) = [764,640]; % RSP
points(2,:) = [848,860]; % V1
points(3,:) = [534,834]; % S1
points(4,:) = [328,748]; % MO
[f1,f2] = butter(2, [2 8]/(35/2), 'bandpass');
for i = 1:4
    trace(i,:) = double(squeeze(Ut1(points(i,1),points(i,2),1:50))'*dV(1:50,:));    
    meanTrace(i,:) = filtfilt(f1,f2,trace(i,:));
end
%%
color1 = cbrewer2('qual','Set1',9);
hs1ab = figure('Renderer', 'painters', 'Position', [100 100 800 300]);
subplot(1,2,1)
scale = 1;
im1 = imagesc(mimgt);
colormap(gray)
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,scale,'w');
axis off; axis image;
for i = 1:4
    hold on;
    scatter(points(i,2),points(i,1),24,color1(i+1,:),'filled');
end

subplot(1,2,2)
plot(t,wheelE/10000+5000,'lineWidth',1,'color','k');
hold on;
plot(t,licking*20+4000,'lineWidth',1,'color',color1(1,:));
hold on;
plot([lick_t,lick_t],[0,1]*400+4000,'k');
for i = 1:4
    hold on;
    plot(t,trace(i,:)*90000+1000*(i-1),'lineWidth',1,'color',color1(i+1,:));
end
hold on;
plot([484,484],[1000, 1900],'r');
yticks([0:1000:5000])
yticklabels({'RSP','VISp','SSp','MO','licking','wheel'})
xlim([476,486]);
xticks([476:2:486])
xlabel('Time (s)');
%%
print(hs1ab, fullfile(save_folder,'FigS1ab_example_trace.pdf'),...
    '-dpdf', '-bestfit', '-painters');
