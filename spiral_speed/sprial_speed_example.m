githubdir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir, 'spikes')))
addpath(genpath(fullfile(githubdir, 'Pipelines')))
addpath(genpath(fullfile(githubdir, 'widefield')))
addpath(genpath(fullfile(githubdir, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
%% SVD, plot example trace, overlay with pupil and stim time
% specify folder for widefield image data
mn = 'ZYE_0012';
td = '2020-10-16';
en = 5;
serverRoot = expPath(mn, td, en);
wf_svd;
%%
figure;
imagesc(mimg)
axis image;
%%
spiral_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals\ZYE_0012\2020-10-16\5';
load(fullfile(spiral_folder,'ZYE_0012_20201016_5_filtered_spirals.mat'));
%%
grid_x = [250,350];
grid_y = [350,550];
index1 = (filteredSpirals2(:,2)>grid_x(1) & filteredSpirals2(:,2)<grid_x(2)...
    &filteredSpirals2(:,1)>grid_y(1) & filteredSpirals2(:,1)<grid_y(2)...
    &filteredSpirals2(:,3)>50);
%%
filteredSpirals3 = filteredSpirals2(index1,:);
%%
lowpass = 0; % 1 is low pass filter <0.2Hz, 0 is bandpass filter between 2-8Hz
Fs = 35; % frame sampling rate
scale = 8;
%%
% [trace2d,traceAmp,tracePhase] = spiralPhaseMap(U1,dV1,t1,lowpass);
[~,~,tracePhase_raw] = spiralPhaseMap(U(1:scale:end,1:scale:end,1:50),dV(1:50,:),t,lowpass);
%%
% padding = 30;
% halfpadding = padding/2;
tracePhase_raw = permute(tracePhase_raw,[2,3,1]);
% [tracePhase_padded] = padZeros(tracePhase_raw,halfpadding);
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 1100 800]);
% for iframe = 1:size(filteredSpirals3,1)
iframe = 2603;
pixSize = 3.45/1000/0.6*3;  % mm / pix
color1 = cbrewer2('seq','Reds',15);
frame = filteredSpirals3(iframe,5);
spiral_radius = (filteredSpirals3(iframe,3))/scale;
center = round([filteredSpirals3(iframe,2)/scale,filteredSpirals3(iframe,1)/scale]);

frame1 = squeeze(tracePhase_raw(:,:,frame));
frame2 = squeeze(tracePhase_raw(:,:,frame+1));

phase_circle1 = []; phase_circle2 = [];
phase_circle1a = []; phase_circle2a = [];
phase_circle1b = []; phase_circle2b = [];
angle_offset = []; distance_offset = [];
angle = pi/6:pi/6:2*pi; % sample evenly around a circle
count2 = 1;
for radius = 1:1:spiral_radius
    for k = 1:numel(angle)
        angle1 = angle(k);
        x(k,count2) = center(2)+(radius*cos(angle1));
        y(k,count2) = center(1)+(radius*sin(angle1));
        phase_circle1(k,count2) = frame1(round(y(k,count2)),round(x(k,count2)));
        phase_circle2(k,count2) = frame2(round(y(k,count2)),round(x(k,count2)));
    end
    count2 = count2+1;
end

angle_offset = [];
% phase angle
for kk = 1:size(phase_circle1,2)
    phase_circle1a = phase_circle1(:,kk);
    phase_circle2a = phase_circle2(:,kk);
    phase_circle1b(:,kk) = wrapAngle(phase_circle1a);
    phase_circle2b(:,kk) = wrapAngle(phase_circle2a);
end
phase_circle1b(abs(phase_circle1b(:))<0.0001) = nan;
phase_circle1b(abs(phase_circle1b(:)-2*pi)<0.0001) = nan;
phase_circle2b(abs(phase_circle2b(:))<0.0001) = nan;
phase_circle2b(abs(phase_circle2b(:)-2*pi)<0.0001) = nan;
angle_offset = phase_circle2b-phase_circle1b;
angle_offset = wrapTo2Pi(angle_offset);
angle_offset(angle_offset(:)>pi) = angle_offset(angle_offset(:)>pi)-pi;
mean_angle_offset = mean(angle_offset,1);
angular_velocity = angle_offset*35;
mean_angular_velocity = mean_angle_offset*35;
%
radius_all = 1:spiral_radius;
radius_all1 = radius_all*pixSize*scale;
kk = numel(radius_all1);

ax1 = subplot(2,4,1);
imagesc(frame1);
colormap(ax1,colorcet('C06'));
axis image;
axis off;
hold on;
scatter(center(2),center(1),12,'w','filled');
hold on;
scatter(x(:,kk),y(:,kk),24,'^','MarkerFaceColor','k','MarkerEdgeColor','none');
for k = 1:numel(angle)
    angle1 = angle(k);
    x1 = center(2)+((spiral_radius+2)*cos(angle1));
    y1 = center(1)+((spiral_radius+2)*sin(angle1));
    text(x1,y1,num2str(k),'fontsize',10);
end
title(['frame' num2str(frame)]);

% scale bar
n1 = round(2/(scale*pixSize));
hold on;
plot([10,10+n1],[60,60],'k','lineWidth',2);

ax2 = subplot(2,4,2);
imagesc(frame2);
colormap(ax2,colorcet('C06'));
axis image;
axis off;
hold on;
scatter(center(2),center(1),12,'w','filled');
hold on;
scatter(x(:,kk),y(:,kk),24,'x','MarkerFaceColor','k','MarkerEdgeColor','k');
for k = 1:numel(angle)
    angle1 = angle(k);
    x1 = center(2)+((spiral_radius+2)*cos(angle1));
    y1 = center(1)+((spiral_radius+2)*sin(angle1));
    text(x1,y1,num2str(k),'fontsize',10);
end
title(['frame' num2str(frame+1)]);

ax3 = subplot(2,4,3);
scatter(1:numel(angle),phase_circle1b(:,kk),24,'^','MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor','none');
hold on;
scatter(1:numel(angle),phase_circle2b(:,kk),24,'x','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on;        
diff1 = phase_circle2b(:,kk)-phase_circle1b(:,kk);
% quiver(1:numel(angle),phase_circle1b(:,kk)',zeros(size(angle)),diff1','off','r')
a = [1:numel(angle)]';
b = [phase_circle1b(:,kk)];c = [phase_circle2b(:,kk)];
hold on; plot([a,a]',[b,c]','r');
ylabel('phase angle wrapped');
xticks([1 2 3 4 5 6 7 8 9 10 11 12]);
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12'});
ylim([-2.5,6]);
xlabel('circular sampling points');

ax4 = subplot(2,4,4);      
diff1 = phase_circle2b(:,kk)-phase_circle1b(:,kk);
scatter(1:numel(angle),diff1,12,'MarkerFaceColor','k','MarkerEdgeColor','none');
ylabel('Angular diff');
xticks([1 2 3 4 5 6 7 8 9 10 11 12])
xticklabels({'1','2','3','4','5','6','7','8','9','10','11','12'})
ylim([15/35 40/35]);
xlabel('circular sampling points');
%
ax7 = subplot(2,4,5);
hold off;
imagesc(frame1);
colormap(ax7,colorcet('C06'));
axis image;
axis off;
hold on;
scatter(center(2),center(1),12,'w','filled');
angle2 = 0:pi/36:2*pi;
radius_all = 1:1:spiral_radius;
cgreys = cbrewer2('seq','greys',numel(radius_all)+4);
count2 = 1;
for radius = 1:1:spiral_radius
    for k = 1:numel(angle2)
        angle3 = angle2(k);
        x1(k,count2) = center(2)+(radius*cos(angle3));
        y1(k,count2) = center(1)+(radius*sin(angle3));
    end
    count2 = count2+1;
end
for radius = 1:1:spiral_radius
    hold on
    plot(x1(:,radius),y1(:,radius),'color',cgreys(radius+4,:),'lineWidth',0.5);
end
hold on;
grid_x = [200,350];
grid_y = [362,512];
v = [362,200;512,200;512,350;362,350]/scale;
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,...
    'EdgeColor','k','FaceColor','none','LineWidth',1);
hold on; 
title(['frame' num2str(frame)]);
%
ax72 = subplot(2,4,6);
imagesc(frame1);
colormap(ax72,colorcet('C06'));
axis image;
axis off;
hold on;
scatter(center(2),center(1),12,'w','filled');
angle2 = 0:pi/36:2*pi;
radius_all = 1:1:spiral_radius;
cgreys = cbrewer2('seq','greys',numel(radius_all)+4);
count2 = 1;
for radius = 1:1:spiral_radius
    for k = 1:numel(angle2)
        angle3 = angle2(k);
        x1(k,count2) = center(2)+(radius*cos(angle3));
        y1(k,count2) = center(1)+(radius*sin(angle3));
    end
    count2 = count2+1;
end
for radius = 1:1:spiral_radius
    hold on
    plot(x1(:,radius),y1(:,radius),'color',cgreys(radius+4,:),'lineWidth',0.5);
    hold on;
    scatter(x(:,radius),y(:,radius),24,'MarkerFaceColor',cgreys(radius+4,:),'MarkerEdgeColor','none');
end
ylim([200/scale,350/scale]);
xlim([362/scale,512/scale]);
% scale bar
n2 = round(0.5/(scale*pixSize));
hold on;
plot([362/scale+1,362/scale+1+n2],[350/scale-1,350/scale-1],'k','lineWidth',2);
%
ax8 = subplot(2,4,7);
hold off;
for radiusi = 1:numel(radius_all1)
    plot(radius_all1(radiusi),angular_velocity(:,radiusi),'-o','color',cgreys(radiusi+4,:),...
    'MarkerSize',3,'MarkerFaceColor',cgreys(radiusi+4,:),'MarkerEdgeColor','none');
    % scatter(raidus_all,angle_offset(angle1,:),'MarkerFaceColor',color1(angle1+2,:),'MarkerEdgeColor','none');
    hold on;
end
hold on;
plot(radius_all1,mean_angular_velocity,'-o','color',[1,0,0],...
'MarkerSize',3,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none');
xlabel('radius (mm)');
ylabel('angluar velocity(rad/s)');
xlim([0,1.6]);
xticks([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6]);
xticklabels({'0' '0.2' '0.4' '0.6' '0.8' '1.0' '1.2' '1.4' '1.6'});


ax9 = subplot(2,4,8);
for kk = 1:numel(radius_all)
    distance_offset(:,kk) = (angle_offset(:,kk)./(2*pi)*(2*pi*radius_all(kk)*pixSize*scale))*35;
end
for radiusi = 1:numel(radius_all1)
    plot(radius_all1(radiusi),distance_offset(:,radiusi),'-o','color',cgreys(radiusi+4,:),...
        'MarkerSize',3,'MarkerFaceColor',cgreys(radiusi+4,:),'MarkerEdgeColor','none');
    % scatter(raidus_all,angle_offset(angle1,:),'MarkerFaceColor',color1(angle1+2,:),'MarkerEdgeColor','none');
    hold on;
end
mean_distance_offset = mean(distance_offset,1);
plot(radius_all1,mean_distance_offset,'-o','color',[1,0,0],...
'MarkerSize',3,'MarkerFaceColor',[1,0,0],'MarkerEdgeColor','none');
xlabel('radius (mm)');
ylabel('velocity (mm/s)');
xlim([0,1.6]);
xticks([0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6]);
xticklabels({'0' '0.2' '0.4' '0.6' '0.8' '1.0' '1.2' '1.4' '1.6'});
%%
print(h1,['spiral_speed_frame' num2str(frame) '-4'], '-dpdf', '-bestfit', '-painters');
close all;
%%
function ph2= wrapAngle(ph)
phdiff = angdiff(ph);
ph2(1) = ph(1);
for i = 2:numel(ph)                
    ph2(i) = [ph2(i-1)+phdiff(i-1)];
end 
end