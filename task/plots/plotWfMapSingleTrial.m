function h1 = plotWfMapSingleTrial(data_folder,T1,session,trial_indx)
%% load atlas brain horizontal projection and outline 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
downscale = 16;
BW3 = BW(1:downscale:end,1:downscale:end);
%% session info
mn = char(T1.MouseID{session});
tda = T1.date(session);  
en = T1.folder(session); 
block_en = T1.block_folder(session); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%%
% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
expDir = dir(fullfile(session_root,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
% get photodiode time
win = [0,3000];
allPD2 = getPhotodiodeTime(session_root,win);
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1)';
allPD2 = allPD2(1:ntrial);
allPD2 = allPD2(norepeat_indx);
%% register to atlas
load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(dV(1:50,:));
%% get time index for all events
for i = 1:numel(allPD2)
    ta  = t-allPD2(i);
    indx2(i,1) = find(ta>0, 1,'first');
    t2(i,1) = t(indx2(i,1));
end
%% plotting params
lineColor = 'k';
hemi = [];
scale3 = 5/16;
lineColor = 'k'; lineColor1 = 'w';
tstep = 0:1/35:18/35;
tstep = round(tstep*1000);
cmax = 0.03; cmin = -0.03;
%% plot
h1 = figure('Renderer', 'painters', 'Position', [50 50 950 300]);
indx_t = indx2(trial_indx);                                                      % stimOn time index in time stamp
Ut1a = double(reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3)));
trace = Ut1a*dV1(:,indx_t:indx_t+32);
trace2d = reshape(trace,size(Ut1,1),size(Ut1,2),size(trace,2))./mimgt;
for i = 1:18
    ax4 = subplottight(1,19,i);
    im1 = imagesc(squeeze(trace2d(:,:,i)));
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
    colormap(ax4,'parula');
    caxis([cmin,cmax]);
    title([num2str(tstep(i)) 'ms']);
end    
cb4 = subplottight(1,19,19);
im1 = imagesc(squeeze(trace2d(:,:,19)),'visible','off');
colormap(cb4,'parula');
caxis([cmin,cmax]);
axis off;
cb = colorbar;
a =  cb.Position; %gets the positon and size of the color bar
set(cb,'Position',[a(1) a(2)+0.45 a(3) a(4)/8])% To change size