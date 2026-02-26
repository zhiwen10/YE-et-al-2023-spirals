function [h1,meanTrace2,tracePhase] = plotMeanTracePhase4(wf_mean,BW2,maskPath,st,atlas1)
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/2;
% v = [60 90; 100 90;100 130;60 130];
v = [80,50;130 50;130 100;80 100]*4;
f = [1,2,3,4];

wf1 = reshape(wf_mean,size(wf_mean,1)*size(wf_mean,2),size(wf_mean,3));
meanTrace = wf1 -mean(wf1 ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
meanTrace2 = meanTrace';
meanTrace2 = reshape(meanTrace2,size(wf_mean,1),size(wf_mean,2),size(wf_mean,3));

traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = tracePhase';
tracePhase = reshape(tracePhase, size(wf_mean,1),size(wf_mean,2),size(wf_mean,3));
mina = min(wf_mean(:));
maxa = max(wf_mean(:));
%%
% figure;
% ax2 = subplottight(1,1,1);
% im_phase = imagesc(tracePhase(:,:,76));
% colormap(ax2,colorcet('C06'));
% axis image;
% set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
% hold on;
% plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
% th2 = 1:5:360; 
% center = [77,104];
% px1 = center(2);
% py1 = center(1);
% r = 20;
% cx2 = (r*cosd(th2)+px1);
% cy2 = (r*sind(th2)+py1);
% cx3 = round(cx2); cy3 = round(cy2);
% hold on;
% color1 = 'w';
% plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
% hold on;
% scatter(center(2),center(1),8,color1,'filled');
%%
h1 = figure('Renderer', 'painters', 'Position', [50 50 1200 400]);
for i = 1:12
    ax1 = subplottight(4,12,i);
    im_raw = imagesc(squeeze(wf_mean(:,:,70+i)));
    % im_raw = imagesc(squeeze(meanTrace2(:,:,66+i)));          
    % im_raw = imagesc(squeeze(wf1(:,:,66+i)));
    caxis([0.8*mina,0.8*maxa]);
    colormap(ax1,parula);
    axis off; 
    axis image;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    pos = get(ax1,'position');
    pos(3) = pos(3)*0.95;
    pos(4) = pos(4)*0.95;
    set(ax1,'position',pos);    

    ax2 = subplottight(4,12,12+i);
    im_phase = imagesc(tracePhase(:,:,70+i));
    colormap(ax2,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    hold on;
    patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);
    caxis([-pi,pi]);
    %
    ax3(i) = subplottight(4,12,12*2+i);
    ax3(i).Position(1) = ax3(i).Position(1);
    ax3(i).Position(2) = ax3(i).Position(2);
    ax3(i).Position(3) = ax3(i).Position(3)-0.02;
    ax3(i).Position(4) = ax3(i).Position(4)-0.02;
    im_raw1 = imagesc(squeeze(wf_mean(:,:,70+i)));
    colormap(ax3(i),parula);  
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    caxis([0.8*mina,0.8*maxa]);
    set(im_raw1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    xlim([80,130]*4);
    ylim([50,100]*4);
    
    ax4(i) = subplottight(4,12,12*3+i);
    ax4(i).Position(1) = ax4(i).Position(1);
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.02;
    ax4(i).Position(4) = ax4(i).Position(4)-0.02;
    im1 = imagesc(tracePhase(:,:,70+i));
    colormap(ax4(i),colorcet('C06'));  
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    axis image; axis off;
    caxis([-pi,pi]);
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    xlim([80,130]*4);
    ylim([50,100]*4);
end
%%
ax1 = subplottight(3,12,12);
im_raw = imagesc(squeeze(wf_mean(:,:,70+12)));
caxis([0.8*mina,0.8*maxa]);
axis off; 
axis image;
set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
cb = colorbar;
cb.Position(4) = cb.Position(4)/2;
pos = get(ax1,'position');
pos(3) = pos(3)*0.95;
pos(4) = pos(4)*0.95;
set(ax1,'position',pos); 