function [count1]= plot_single_pixel2(ops,sp,mean_var_t, mimgt,mimg, BW1, point, point_t,...
    data_predict, data_raw,WF2ephysT,epochs,Axy_all,frame,filename)
%%
scale = 4;
h1 = figure('Renderer', 'painters', 'Position', [100 100 1000 600]);
ax1 = subplot(4,4,1);
im1 = imagesc(mimgt);
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
colormap(ax1,gray);
axis image; axis off;
mxRange = prctile(mimgt(:), 99.5);
caxis([0,mxRange])
colorbar;
hold on;
scatter(point_t(:,2),point_t(:,1),6,'r','filled')
%%
color1 = cbrewer2('qual','Set2',8);
% cmap1 = flipud(cbrewer2('div', 'RdYlBu', 32,'linear'));
ax2 = subplot(4,4,5);
im2 = imagesc(mean_var_t);
set(im2, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
% caxis([0,1]);
axis image; axis off;
cmap1 = flipud(inferno);
colormap(ax2,cmap1);
colorbar;

point1 = round(point/scale);
mimg_scale = mimg(1:scale:end,1:scale:end);
subplot(4,4,[2:4])
trace_predict = squeeze(data_predict.trace2d1(:,point1(1,1),point1(1,2)))./mimg_scale(point1(1),point1(2));
trace_raw = squeeze(data_raw.trace2d1(:,point1(1,1),point1(1,2)))./mimg_scale(point1(1),point1(2));
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_raw*100,'lineWidth',1,'color',color1(5,:)); % GREEN
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_predict*100,'lineWidth',1,'color',color1(8,:)); % GRAY
xlim([WF2ephysT(epochs(1)),WF2ephysT(epochs(end))]);
xlabel('Time (s)');
ylabel('df/f (%)');
hold on;
xline(WF2ephysT(epochs(frame(1))),'--');
xline(WF2ephysT(epochs(frame(end))),'--');
hold off;
%
t_sample = (epochs(end)-epochs(1));
t_length = t_sample/35;

subplot(4,4,[6:8])
trace_predict = squeeze(data_predict.trace2d1(:,point1(2,1),point1(2,2)))./mimg_scale(point1(1),point1(2));
trace_raw = squeeze(data_raw.trace2d1(:,point1(2,1),point1(2,2)))./mimg_scale(point1(1),point1(2));
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_raw*100,'lineWidth',1,'color',color1(5,:)); % GREEN
hold on;
plot(WF2ephysT(epochs(1:end-1)),trace_predict*100,'lineWidth',1,'color',color1(8,:)); % GRAY
xlim([WF2ephysT(epochs(1)),WF2ephysT(epochs(end))]);
xlabel('Time (s)');
ylabel('df/f (%)');
hold on;
xline(WF2ephysT(epochs(frame(1))),'--');
xline(WF2ephysT(epochs(frame(end))),'--');
hold off;
%
t_sample = (epochs(end)-epochs(1));
t_length = t_sample/35;

subplot(4,4,[10:12])
% ishank = 1:4;
% incl = (sp.spikeAmps>20 & ismember(sp.clu,sp.gcluster) & ismember(sp.spikeSites,ops.chanMap1(:,ishank)));
color1 = colorcet('C06','N',8);
count1 = 0;
for i = 1:numel(sp.gcluster)
    clear incl cluster_t cluster_depth
%     incl = (sp.clu==sp.gcluster(i) & sp.spikeTimes>WF2ephysT(epochs(1)) ...
%         & sp.spikeTimes<WF2ephysT(epochs(end))...
%         & ismember(sp.spikeSites,ops.chanMap1(:,ishank)));
    incl = (sp.clu==sp.gcluster(i) & sp.spikeTimes>WF2ephysT(epochs(1)) ...
    & sp.spikeTimes<WF2ephysT(epochs(end)));
    cluster_t = sp.spikeTimes(incl);
    cluster_depth = sp.spikeDepths(incl);
    mean_rate = numel(cluster_t)./t_length;
    if mean_rate<20
        randi = randperm(8,1);
        % randi = mod(i,7)+1;
        plot([cluster_t,cluster_t]',[cluster_depth,cluster_depth+100]','Color',color1(randi,:),'lineWidth',1);
        hold on;
        count1 = count1+1;
    end
end
% ylim([-200,1600]);
xlim([WF2ephysT(epochs(1)),WF2ephysT(epochs(end))]);

%%
BW2 = BW1(1:8:end,1:8:end);
mean_var_t2 = mean_var_t(1:8:end,1:8:end);
indx = find(mean_var_t2<=0.4 & not(BW2));
ax10 = subplot(4,4,14:16);
theta = angle(Axy_all);
theta(:,indx) = [];
[S s] = circ_var(theta, [], [], 2);
c = 1-S;
scatter(1:t_sample,c,[],'k','filled','MarkerFaceAlpha',0.4)
xticks([0:17.5:t_sample]);
xticksN = num2cell(0:0.5:t_length);
xticksStr = cellfun(@num2str,xticksN,'UniformOutput',false);
xticklabels(xticksStr);
xlabel('Time (s)');
xlim([1,t_sample]);
ylim([0,1])
%%
print(h1, filename, '-dpdf', '-bestfit', '-painters');
end