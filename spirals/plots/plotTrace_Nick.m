%% download circular color map
githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';          % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'colorcet')));                         % https://colorcet.com/download/index.html
%% data prep only
% qt2 = qta(frame0:last_frame1+10);
% trace_raw2 = trace_raw(:,frame0:last_frame1+10);
% trace_filt2 = trace_filt(:,frame0:last_frame1+10);
% trace_phase2 = trace_phase(:,frame0:last_frame1+10);
% save("epoch_data.mat",'qt2', 'pixel', 'trace_raw2', 'trace_phase2', 'trace_filt2');
%% plotting
% data was upsampled 10 times, so that circular color points are smooth
% load data: qt2, pixel, trace_raw2, trace_phase2, trace_filt2
load("epoch_data.mat");
% define the epoch window (31:end-10)
first_frame2 = 31; last_frame2 = length(qt2)-10;
h1ac = figure('Renderer', 'painters', 'Position', [100 100 200 600]);
ax1 = subplot(1,1,1);
yscale = 150;
for i = 1:size(pixel,1)
    plot(ax1,qt2,trace_raw2(i,:)*yscale-3*(i-1),'k','lineWidth',1);
    hold on; 
    %% convert phase range of[-pi,pi] to [-0.5,0.5] then to [0,1], so it can
    % be interpolated
    color1 = trace_phase2(i,:)/(2*pi)+0.5; 
    % define 100 points in the circular color map
    color_hsv = colormap(ax1,colorcet('C06','N',100));
    num1 = linspace(0,1,100);
    % interpolate the color of phase trace for each time point
    color_hsv1 = interp1(num1,color_hsv,color1);
    qt1a = qt2(1:end-1);qt2a = qt2(2:end);
    % plot circular color by draw lines between 2 points
    trace1 = trace_filt2(i,1:end-1);
    trace2 = trace_filt2(i,2:end);
    for j = 1:numel(qt1a)
    plot(ax1,[qt1a(j),qt2a(j)],[trace1(j)*yscale-3*(i-1),trace2(j)*yscale-3*(i-1)],...
        'lineWidth',2,'color',color_hsv1(j,:));
    hold on; 
    end
    %% find and plot peak points
    trace_filt_temp = trace_filt2(i,:);
    [pks,locs] = findpeaks(trace_filt_temp,'MinPeakProminence',0.001);
    scatter(qt2(locs),trace_filt_temp(locs)*yscale-3*(i-1)+0.2,12,'k','filled');    
    hold on; 
end
ax1.FontSize = 8; 
set(ax1,'YTickLabel',[]);
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(qt2(first_frame2),'--');
line_ax2 = xline(qt2(last_frame2),'--');
