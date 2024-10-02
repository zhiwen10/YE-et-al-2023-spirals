function h1e = plotSpiralDensityAllSessions_freq(freq,data_folder,save_folder)
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%% load total frame numbers from all 15 sessions
load(fullfile(data_folder,'spirals\spirals_density',...
    'histogram_40pixels.mat'),'frame_total');
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
load(fullfile(data_folder,'spirals\spirals_freq\spirals_density',...
    ['histogram_' freq_folder '.mat']));
%%
h1e = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
ax1 = subplot(1,1,1);
hist_bin = 40;
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_total*35;                 % spirals/(mm^2*s)
scatter(ax1,unique_spirals(:,1),unique_spirals(:,2),...
    3,unique_spirals_unit,'filled');
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
colormap(ax1,hot);
% axis off; 
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
cb_size1 =  get(cb1,'Position');                                           % gets the positon and size of the color bar
set(cb1,'Position',...
    [cb_size1(1)+0.07  cb_size1(2) cb_size1(3) cb_size1(4)/2]);            % To change size
title(freq_folder,'interpreter','none');
print(h1e, fullfile(save_folder,['histogram_' freq_folder '.pdf']), ...
    '-dpdf', '-bestfit', '-painters');