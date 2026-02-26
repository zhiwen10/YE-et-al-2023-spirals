function [h1,meanTrace2,tracePhase] = plotMeanTracePhase7(wf_mean2,BW2,maskPath,st,atlas1,spirals)
%%
maskPath{12} = '/997/8/567/688/695/315/453/322/329/'; % SSp-bfd
maskPath{13} ='/997/8/567/688/695/315/453/322/353/'; %SSP-n
maskPath{14} ='/997/8/567/688/695/315/453/322/337/'; %SSp-ll
maskPath{15} ='/997/8/567/688/695/315/453/322/345/'; %SSp-m
maskPath{16} ='/997/8/567/688/695/315/453/322/369/'; %SSp-ul
maskPath{17} ='/997/8/567/688/695/315/453/322/361/'; %SSp-tr
maskPath{18} ='/997/8/567/688/695/315/453/322/182305689/'; %SSp-un
%%
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/2;
% v = [60 90; 100 90;100 130;60 130];
v = [80,50;130 50;130 100;80 100]*4;
f = [1,2,3,4];

wf1 = reshape(wf_mean2,size(wf_mean2,1)*size(wf_mean2,2),size(wf_mean2,3));
meanTrace = wf1 -mean(wf1 ,2);
meanTrace = double(meanTrace)';
[f1,f2] = butter(2,[2,8]/(35/2), 'bandpass');
meanTrace = filtfilt(f1,f2,meanTrace);
meanTrace2 = meanTrace';
meanTrace2 = reshape(meanTrace2,size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));

traceHilbert =hilbert(meanTrace);
tracePhase = angle(traceHilbert);
tracePhase = tracePhase';
tracePhase = reshape(tracePhase, size(wf_mean2,1),size(wf_mean2,2),size(wf_mean2,3));
% mina = min(wf_mean2(:));
% maxa = max(wf_mean2(:));
mina = -0.01; maxa = 0.01;
th2 = 1:5:360; 
%% 
subN = 14;
h1 = figure('Renderer', 'painters', 'Position', [50 50 950 400]);
for i = 1:subN
    spiral_temp = spirals(spirals(:,5) == 68+i,:);
        
    ax1 = subplottight(4,15,i);
    im_raw = imagesc(squeeze(wf_mean2(:,:,68+i)));
    % im_raw = imagesc(squeeze(meanTrace2(:,:,66+i)));          
    % im_raw = imagesc(squeeze(wf1(:,:,66+i)));
    % caxis([0.8*mina,0.8*maxa]);
    caxis([mina,maxa]);
    colormap(ax1,parula);
    axis off; 
    axis image;
    set(im_raw, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    hold on;
    if i==3
        for kk = 12:18
            plotOutline(maskPath(kk),st,atlas1,hemi,scale3,lineColor);
        end
    end
    hold on;
    patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);
    pos = get(ax1,'position');
    pos(3) = pos(3)*0.95;
    pos(4) = pos(4)*0.95;
    set(ax1,'position',pos);    

    ax2 = subplottight(4,15,15+i);
    im_phase = imagesc(tracePhase(:,:,68+i));
    colormap(ax2,colorcet('C06'));
    axis image; axis off;
    set(im_phase, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    hold on;
    if i==3
        for kk = 12:18
            plotOutline(maskPath(kk),st,atlas1,hemi,scale3,lineColor);
        end
    end
    hold on;
    patch('Faces',f,'Vertices',v,'EdgeColor','k','FaceColor','none','LineWidth',0.8);
    caxis([-pi,pi]);
    %
    ax3(i) = subplottight(4,15,15*2+i);
    ax3(i).Position(1) = ax3(i).Position(1)+0.01;
    ax3(i).Position(2) = ax3(i).Position(2);
    ax3(i).Position(3) = ax3(i).Position(3)-0.02;
    ax3(i).Position(4) = ax3(i).Position(4)-0.02;
    im_raw1 = imagesc(squeeze(wf_mean2(:,:,68+i)));
    colormap(ax3(i),parula);  
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if i==3
        for kk = 12:18
            plotOutline(maskPath(kk),st,atlas1,hemi,scale3,lineColor);
        end
    end
    hold on;
    axis image; axis off;
    % caxis([0.8*mina,0.8*maxa]);
    caxis([mina,maxa]);
    set(im_raw1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end 
    xlim([80,130]*4);
    ylim([50,100]*4);
    
    ax4(i) = subplottight(4,15,15*3+i);
    ax4(i).Position(1) = ax4(i).Position(1)+0.01;
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.02;
    ax4(i).Position(4) = ax4(i).Position(4)-0.02;
    im1 = imagesc(tracePhase(:,:,68+i));
    colormap(ax4(i),colorcet('C06'));  
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
    if i==3
        for kk = 12:18
            plotOutline(maskPath(kk),st,atlas1,hemi,scale3,lineColor);
        end
    end
    hold on;
    axis image; axis off;
    caxis([-pi,pi]);
    set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled'); 
    if not(isempty(spiral_temp))
        for kk = 1:size(spiral_temp,1)
            px1 = spiral_temp(kk,1);
            py1 = spiral_temp(kk,2);
            r = spiral_temp(kk,3);
            cx2 = round(r*cosd(th2)+px1);
            cy2 = round(r*sind(th2)+py1);
            hold on;
            if spiral_temp(kk,4) == 1                                                    % counterclockwise, then color white
                color1 = 'w';
            else
                color1 = 'k';                                                      % clockwise, then color black
            end
            plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);          % draw the circle at max radius
            hold on;
            scatter(spiral_temp(kk,1),spiral_temp(kk,2),8,color1,'filled');
        end
    end
    xlim([80,130]*4);
    ylim([50,100]*4);
end
%
ax1 = subplottight(3,15,15);
im_raw = imagesc(squeeze(wf_mean2(:,:,68+12)));
% caxis([0.8*mina,0.8*maxa]);
caxis([mina,maxa]);
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