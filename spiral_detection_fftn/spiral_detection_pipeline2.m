githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')));
addpath(genpath(fullfile(githubdir2, 'Pipelines')));
addpath(genpath(fullfile(githubdir2, 'widefield')));
addpath(genpath(fullfile(githubdir2, 'npy-matlab')));
addpath('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection');
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
% addpath(genpath(fullfile(githubdir2, 'NeuroPattToolbox'))) %https://github.com/BrainDynamicsUSYD/NeuroPattToolbox
% addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
T = readtable('spiralSessions.xlsx');
% codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow1';
% T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
%%
kk = 8;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%% SVD, plot example trace, overlay with pupil and stim time
serverRoot = expPath(mn, td, en);
[U,V,t,mimg] = get_wf_svd1(serverRoot);
dV = [zeros(size(V,1),1) diff(V,[],2)];
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
params.epochL = 100;
params.nt = numel(t);
frameN1 = round((params.nt-70)/100)*100;
params.frameRange = 36:params.epochL:(36+frameN1-1);
params.dThreshold = 15;
params.rsRCheck = 10:10:100;
%% roi grid
% only look around roi with gridsize radius
[params.pgridx,params.pgridy] = meshgrid(21-10:21+10, 21-10:21+10);
%%
U1 = U(1:params.downscale:end,1:params.downscale:end,1:50);
params.xsize = size(U1,1);
params.ysize = size(U1,2);
xsizePadded = params.xsize+params.padding; ysizePadded = params.ysize+params.padding;
[xx,yy] = meshgrid(min(params.rs)+1:params.gridsize:xsizePadded-min(params.rs)-1, ...
    min(params.rs)+1:params.gridsize:ysizePadded-min(params.rs)-1);
params.xx = xx; params.yy = yy;
%% apply mask, this helps speed up spiral detection later
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
mimg2 = padZeros(mimg1,params.halfpadding);
%%
% figure;
% imagesc(mimg2);
% axis image;
% roi = drawpolygon;
% save([mn '_' tdb '_' num2str(en) '_roi'], 'roi');
%%
fname1 = [mn '_' tdb '_' num2str(en) '_roi'];
folder1 = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals\ZYE_0012\2020-10-16\5';
load(fullfile(folder1, [fname1 '.mat']));
tf = inROI(roi,xx(:),yy(:));
params.xxRoi = xx(tf); 
params.yyRoi =yy(tf);
rate = 1;
%%
pwAll = [];
tic
pwAll = [];
pwAll1 = []; pwAll2 = []; pwAll3 = []; pwAll4 = []; pwAll5 = [];
tStart = 1681; tEnd = 1684; % find spirals between time tStart:tEnd
frameStart = find(t>tStart,1,'first'); frameEnd = find(t>tEnd,1,'first');
frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
dV1 = dV(1:50,frameTemp);
t1 = t(frameTemp);      
[trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U1,dV1,t1,params,rate);
tracePhase1 = tracePhase1(:,:,1+35:end-35); % reduce 2*35 frames after filter data 
tracePhase = padZeros(tracePhase1,params.halfpadding); % pad tracephase with edge zeros
nframe = size(tracePhase,3);
for frame = 1:nframe
    A = squeeze(tracePhase(:,:,frame));
    [pwAll1] = spiralAlgorithm(A,params);
    frameN = frame+frameStart-1;
    pwAll2 = checkClusterXY(pwAll1,params.dThreshold);
    % double check the mean cluster centers from the last step are still
    % spiral centers
    [pwAll3] = doubleCheckSpiralsAlgorithm(A,pwAll2,params);
    [pwAll4] = spatialRefine(A,pwAll3,params);
    [pwAll5] = spiralRadiusCheck2(A,pwAll4,params);
    if not(isempty(pwAll5))
        pwAll5(:,end+1) = frameN; 
    end
    pwAll = [pwAll;pwAll5];
end
fprintf('Frame %g/%g; time elapsed %g seconds \n', [frameStart,frameN1, toc])
pwAll(:,1:2) = pwAll(:,1:2)-params.halfpadding;
%%
frame2 = 58744:58944; 
rate1 = 1;
dV1 = dV(1:50,frame2);
t1 = t(frame2);
[trace2d2,traceAmp2,tracePhase2] = spiralPhaseMap4(U(:,:,1:50),dV1,t1,params,rate1);
%%
iframe = 21;
pwi = 8;
% iframe = 25;
% pwi = 13;
tracePhase1 = squeeze(tracePhase2(:,:,iframe));
tracePhase = padZeros(tracePhase1,params.halfpadding); % pad tracephase with edge zeros
[pwAll1] = spiralAlgorithm(tracePhase,params);
%%
figure;
subplot(1,2,1)
imagesc(squeeze(tracePhase2(:,:,iframe)));
colormap(colorcet('C06'));
axis image; axis off;
subplot(1,2,2)
imagesc(squeeze(tracePhase2(:,:,iframe+1)));
colormap(colorcet('C06'));
axis image; axis off;
%%

h = figure('Renderer', 'painters', 'Position', [100 100 1100 600]);
ax3 = subplot(2,3,1);
imagesc(tracePhase);
colormap(ax3,colorcet('C06'));
axis image;
hold on;
scatter(params.xxRoi,params.yyRoi,0.5,'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');
hold on;
px= pwAll1(pwi,1);
py = pwAll1(pwi,2);
scatter(px,py,8,'k','filled');
hold on;
v = [500 350; 600 350; 600 450; 500 450];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
hold on; 
plot([400,515],[700,700],'k')


greys  = cbrewer2('seq','Greys',9);
ax3 = subplot(2,3,2);
rs = params.rs;
th = params.th;
spiralRange = params.spiralRange;
imagesc(tracePhase);
colormap(ax3,colorcet('C06'));
axis image;
hold on;
scatter(params.xxRoi,params.yyRoi,3,'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');
hold on;
scatter(px,py,24,'k','filled','markerEdgeColor','None');
for rn = 1:numel(rs)
    r = rs(rn);
    cx = round(r*cosd(th)+px);
    cy = round(r*sind(th)+py);
    hold on;
    scatter(cx,cy,24,'markerFacecolor',greys(2+2*rn,:),'markerEdgeColor','None');
    hold on;
    plot([cx cx(1)],[cy cy(1)],'color',greys(2+2*rn,:),'LineWidth',2);
end
xlim([500,600]);
ylim([350,450]);

ax4 = subplot(2,3,3);
hold off;
for rn = 1:numel(rs)
    r = rs(rn);
    cx = round(r*cosd(th)+px);
    cy = round(r*sind(th)+py);
    ph = tracePhase(sub2ind(size(tracePhase),cy, cx));
    phdiff = angdiff(ph);
    ph2(1) = ph(1);
    for i = 2:numel(ph)                
        ph2(i) = [ph2(i-1)+phdiff(i-1)];
    end 
    ph3 = abs(ph2-ph2(1));                
    [N,edges] = histcounts(ph,spiralRange);
    AngleRange = abs(ph2(end)-ph2(1));
    %%
    scatter(1:10,ph3,24,'MarkerFaceColor',greys(2+2*rn,:),'MarkerEdgeColor','None');
    hold on;
end
xlabel('Sampling points');
ylabel('Cumulative phase angle');

ax5 = subplot(2,3,4);
imagesc(tracePhase);
colormap(ax5,colorcet('C06'));
axis image;
% hold on;
% scatter(params.xxRoi,params.yyRoi,1,'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');
hold on;
scatter(pwAll1(:,1),pwAll1(:,2),8,'k','filled');

ax6 = subplot(2,3,5);
pwAll2 = checkClusterXY(pwAll1,params.dThreshold);
% double check the mean cluster centers from the last step are still
% spiral centers
[pwAll3] = doubleCheckSpiralsAlgorithm(tracePhase,pwAll2,params);
imagesc(tracePhase);
colormap(ax6,colorcet('C06'));
axis image;
% hold on;
% scatter(params.xxRoi,params.yyRoi,1,'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');
hold on;
scatter(pwAll3(:,1),pwAll3(:,2),8,'k','filled');
%
[pwAll4] = spatialRefine(tracePhase,pwAll3,params);
[pwAll5] = spiralRadiusCheck2(tracePhase,pwAll4,params);
ax7 = subplot(2,3,6);
th2 = 1:5:360; 
imagesc(tracePhase);
colormap(ax7,colorcet('C06'));
axis image;
hold on;
scatter(pwAll5(:,1),pwAll5(:,2),8,'k','filled');
hold on;
for i = 1:size(pwAll5,1)
    px1 = pwAll5(i,1);
    py1 = pwAll5(i,2);
    r = pwAll5(i,3);
    cx2 = round(r*cosd(th2)+px1);
    cy2 = round(r*sind(th2)+py1);
%     hold on;
%     scatter(cx,cy,6,'markerFacecolor','k','markerEdgeColor','None');
    hold on;
    if pwAll5(i,4) == 1 % counterclockwise
        color1 = 'w';
    else
        color1 = 'k';
    end
    plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',1);
end
%%
print(h, 'detection-pipeline-1', '-dpdf', '-bestfit', '-painters');
%%
radiusSize = 3.45/1000/0.6*3;  % mm / pix
  % that's 3.45 um/pix (https://www.baslerweb.com/en/products/cameras/area-scan-cameras/ace/aca2440-75um/)
  % converted to mm
  % adjusted for 0.6x magnification
  % adjusted for 3x3 binning
  % (the answer is 17.3 um/pix)
% pixArea = pixSize^2; %mm^2