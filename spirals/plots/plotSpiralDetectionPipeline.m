function hs2 = plotSpiralDetectionPipeline(T,data_folder,save_folder)
%% use ZYE12 as an example
kk = 7;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%% load SVD data
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);
dV = [zeros(size(V,1),1) diff(V,[],2)];
%% set params for spiral detection
freq = [2,8];                                                              % data filtering frequency range
params = setSpiralDetectionParams(U,t);                                    % set detection parameters
fname1 = [mn '_' tdb '_' num2str(en) '_roi'];                              % apply mask, this helps speed up spiral detection later
load(fullfile(data_folder,'spirals','full_roi', [fname1 '.mat']));
tf = inROI(roi,params.xx(:),params.yy(:));
params.xxRoi = params.xx(tf); 
params.yyRoi = params.yy(tf);
%% look at these frame range for example sprials
frame2 = 58744:58944; 
rate1 = 1;
dV1 = dV(1:50,frame2);
t1 = t(frame2);
[trace2d2,traceAmp2,tracePhase2] = ...
    spiralPhaseMap_freq(U(:,:,1:50),dV1,t1,params,freq,rate1);
%% frame 21 has a beautiful sprial as a good example
iframe = 21;
pwi = 8;
tracePhase1 = squeeze(tracePhase2(:,:,iframe));
tracePhase = padZeros(tracePhase1,params.halfpadding);                     % pad tracephase with edge zeros
[pwAll1] = spiralAlgorithm(tracePhase,params);                             % detect sprials on this frame
%%
hs2 = figure('Renderer', 'painters', 'Position', [100 100 1100 600]);
ax1 = subplot(2,3,1);                                                      % plot padded example frame
imagesc(tracePhase);
colormap(ax1,colorcet('C06'));
axis image;
hold on;
scatter(params.xxRoi,params.yyRoi,0.5,...
    'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');             % overlaid with searching grids 
hold on;
px= pwAll1(pwi,1);
py = pwAll1(pwi,2);
scatter(px,py,8,'k','filled');                                             % example spiral center 
hold on;
v = [500 350; 600 350; 600 450; 500 450];                                  % zoom in patch
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
hold on; 
plot([400,515],[700,700],'k')

greys  = cbrewer2('seq','Greys',9);
ax2 = subplot(2,3,2);                                                      % plot phase map within zoom in path
rs = params.rs;
th = params.th;
spiralRange = params.spiralRange;
imagesc(tracePhase);
colormap(ax2,colorcet('C06'));
axis image;
hold on;
scatter(params.xxRoi,params.yyRoi,3,...
    'markerEdgecolor',[0.5,0.5,0.5],'markerFaceColor','None');             % plot searching grids
hold on;
scatter(px,py,24,'k','filled','markerEdgeColor','None');                   % example spiral center
for rn = 1:numel(rs)
    r = rs(rn);
    cx = round(r*cosd(th)+px);
    cy = round(r*sind(th)+py);
    hold on;
    scatter(cx,cy,24,'markerFacecolor',greys(2+2*rn,:),...
        'markerEdgeColor','None');                                         % sampling points on a circle
    hold on;
    plot([cx cx(1)],[cy cy(1)],'color',greys(2+2*rn,:),'LineWidth',2);     % join the points on the circle   
end
xlim([500,600]);
ylim([350,450]);

ax3 = subplot(2,3,3);                                                      % plot cumulative phase angles along a circle
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
    scatter(1:10,ph3,24,'MarkerFaceColor',greys(2+2*rn,:),...
        'MarkerEdgeColor','None');
    hold on;
end
xlabel('Sampling points');
ylabel('Cumulative phase angle');

ax4 = subplot(2,3,4);                                                      % plot all candidate spirals on a single frame 
imagesc(tracePhase);
colormap(ax4,colorcet('C06'));
axis image;
hold on;
scatter(pwAll1(:,1),pwAll1(:,2),8,'k','filled');

ax5 = subplot(2,3,5);
pwAll2 = checkClusterXY(pwAll1,params.dThreshold);                         % cluster nearby duplicte spiral centers
[pwAll3] = doubleCheckSpiralsAlgorithm(tracePhase,pwAll2,params);          % double check if it is still a spiral
imagesc(tracePhase);
colormap(ax5,colorcet('C06'));
axis image;
hold on;
scatter(pwAll3(:,1),pwAll3(:,2),8,'k','filled');
[pwAll4] = spatialRefine(tracePhase,pwAll3,params);                        % refined search of sprial centers
[pwAll5] = spiralRadiusCheck2(tracePhase,pwAll4,params);                   % check spiral radius and direction
    
ax6 = subplot(2,3,6);                                                      % plot final spirals with radius and direction
th2 = 1:5:360; 
imagesc(tracePhase);
colormap(ax6,colorcet('C06'));
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
    hold on;
    if pwAll5(i,4) == 1                                                    % counterclockwise, then color white
        color1 = 'w';
    else
        color1 = 'k';                                                      % clockwise, then color black
    end
    plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',1);          % draw the circle at max radius
end
print(hs2, fullfile(save_folder,'FigS2_detection-pipeline'),...
    '-dpdf', '-bestfit', '-painters');                                     % save as pdf file
