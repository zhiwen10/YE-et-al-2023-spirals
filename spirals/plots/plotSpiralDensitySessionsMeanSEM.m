function hs7cd = plotSpiralDensitySessionsMeanSEM(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
load(fullfile(data_folder,'spirals','spirals_density',...
    'spiralDensityLinePerSession.mat'))
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
load(fullfile(data_folder,'spirals','spirals_density',...
    'histogram_40pixels.mat'));
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
mean_g7 = mean(max_density(1:11));
std_g7 = std(max_density(1:11));
mean_g8 = mean(max_density(12:15));
std_g8 = std(max_density(12:15));
[hh,p,ci,stats] = ttest2(max_density(1:11),max_density(12:15));
%%
session_total = 15;
mean_spirals = mean(count_sample,2);
std_spirals = std(count_sample,[],2)./sqrt(session_total);
%%
hs7cd = figure;
clear unique_spiral_unit
hist_bin = 40;
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_total*35;                 % spirals/(mm^2*s)
ax1 = subplot(1,2,1);
scatter(unique_spirals(:,1),unique_spirals(:,2),3,unique_spirals_unit,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
colormap(ax1,hot);
axis image;
xlim(ax1,[0,1140]);
ylim(ax1,[0,1320]);
xtick = 0:200:1000;
xtick_string = string(xtick);
xticks(xtick)
xticklabels(xtick_string)
ytick = 0:200:1200;
ytick_string = string(ytick);
yticks(ytick)
yticklabels(ytick_string)
cmax = max(unique_spirals_unit);
caxis([0,cmax]);
cb1 = colorbar;
cmax = max(unique_spirals_unit(:));
cb1.Ticks = [0,cmax];
cb1.TickLabels = {'0',num2str(round(cmax*10)/10)};
cb_size1 =  get(cb1,'Position'); %gets the positon and size of the color bar
set(cb1,'Position',[cb_size1(1)+0.07  cb_size1(2) cb_size1(3) cb_size1(4)/2]);% To change size
% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
hold on;
scatter(points(:,1),points(:,2),8,'g','filled');
% interp histgram counts 
% F = scatteredInterpolant(unique_spirals(:,1),unique_spirals(:,2),unique_spirals_unit );
% xx = 1:1140; yy = 1:1320;
% [xxq,yyq] = meshgrid(xx,yy);
% vq1 = F(xxq,yyq);
% for i = 1:size(points,1)
%     count_sample(i,1) = vq1(points(i,2),points(i,1));
% end
subplot(1,2,2);
points1 = points-points(1,:);
points_line = vecnorm(points1,2,2);
shadedErrorBar(points_line(:,1),mean_spirals,std_spirals,'lineProps','g');
xlim([0,800]);
ylim([0,2]);
yticks([0,0.5,1,1.5,2]);
yticklabels({'0' '0.5' '1' '1.5' '2'});
xticks([0 200 400 600 800]);
xticklabels({'0','2','4','6','8'});
ylabel('sprials/(mm^2*s)');  
%%
print(hs7cd, fullfile(save_folder,'Figs7cd_spirals_density_sessions'),...
    '-dpdf', '-bestfit', '-painters');