githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines')))
addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox'))) %https://github.com/BrainDynamicsUSYD/NeuroPattToolbox
addpath(genpath(fullfile(githubdir2, 'widefield')))
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\spiralDetection'))
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
rf_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\receptive_field_mapping';
scale = 1;
session_all = find(T.use);
session_total = numel(session_all);
%%
params.downscale = 1; % downscale factor for each frame to process faster, by sacrificing resolution
params.lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
params.Fs = 35; % frame sampling rate
params.halfpadding = 120; % half padding shoudl be the max spiral detection radius 
params.padding = 2*params.halfpadding;
params.th = 1:36:360; 
params.rs = 10:5:20;
params.gridsize = 10;
params.spiralRange = linspace(-pi,pi,5);
params.gsmooth = 0;
params.epochL = 1000;
params.dThreshold = 15;
params.rsRCheck = 10:10:100;
rate = 1;
%% roi grid
% only look around roi with gridsize radius
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);
%%
kk = 15;
mn = T.MouseID{session_all(kk)};
tda = T.date(session_all(kk));
en = T.folder(session_all(kk));
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
serverRoot = expPath(mn, td, en);
mimg = readNPY(fullfile(serverRoot, 'blue','meanImage.npy'));
fname = [mn '_' tdb '_' num2str(en) '_roi'];
figure; 
subplot(1,2,1);
imagesc(mimg);
axis image;
roi_ml = [1,550];
roi_ap = [200,540];
% roi_ml = [1,512];
% roi_ap = [200,480];
hold on;
rectangle('Position',[roi_ml(1),roi_ap(1),roi_ml(2)-roi_ml(1),roi_ap(2)-roi_ap(1)],'EdgeColor','w');
text(100,100,mn,'color','w','interpreter','None')
subplot(1,2,2); 
imagesc(mimg(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2)));
axis image;
save(fname,'roi_ap','roi_ml');
%%
mimg = mimg(roi_ap(1):roi_ap(2),roi_ml(1):roi_ml(2));
params.xsize = size(mimg,1);
params.ysize = size(mimg,2);
xsizePadded = params.xsize+params.padding; ysizePadded = params.ysize+params.padding;
[xx,yy] = meshgrid(min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
    min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
params.xx = xx; params.yy = yy;
%% apply mask, this helps speed up spiral detection later
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
mimg2 = padZeros(mimg1,params.halfpadding);
mimg2_bw = logical(mimg2);
figure; 
ax1 = imagesc(mimg2_bw);
hold on;
scatter(yy,xx,'r','filled');
%%
tf = zeros(size(xx));
for i = 1:numel(xx(:))
    tf(i) = mimg2_bw(xx(i),yy(i));
end
tf = logical(tf);
%%
params.xxRoi = yy(tf); 
params.yyRoi =xx(tf);
