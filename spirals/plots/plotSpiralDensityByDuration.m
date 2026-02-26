function hs6cd = plotSpiralDensityByDuration(data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
%%
load(fullfile(data_folder,'spirals','spirals_density_duration',...
    'frame_total.mat'));
duration = 1:7;
duration_t = round(duration/35*1000);
hist_bin = 40;
count1 = 1;
hs6cd = figure('Renderer', 'painters', 'Position', [100 100 1300 700]);
for k  = 1:7
    clear unique_spiral_unit
    duration_i = duration(k);
    load(fullfile(data_folder,'spirals','spirals_density_duration',...
        ['histogram_' num2str(duration_i) 'frame.mat']));                  % load histogtam counts
    unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./frame_total*35;             % spirals/(mm^2*s)
    ax1 = subplot(2,7,count1);
    scatter(unique_spirals(:,1),unique_spirals(:,2),...
        3,unique_spirals_unit,'filled');
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
    xticklabels({})
    ytick = 0:200:1200;
    yticks(ytick)
    yticklabels({})
    cb1 = colorbar;
    cmax = max(unique_spirals_unit(:));
    cb1.Ticks = [0,cmax];
    cb1.TickLabels = {'0',num2str(round(cmax*10)/10)};
    cb_size1 =  get(cb1,'Position'); %gets the positon and size of the color bar
    set(cb1,'Position',[cb_size1(1)+0.05  cb_size1(2) cb_size1(3) cb_size1(4)/2]);% To change size
    title([num2str(duration_t(k)) 'ms']);
    %% draw a line 
    x = [105,520]; y = [780,225];
    xq = 105:520;
    vq = interp1(x,y,xq);
    points = [xq',round(vq)'];
    hold on;
    scatter(points(:,1),points(:,2),5,'g','filled');
    %% interp histgram counts 
    F = scatteredInterpolant(unique_spirals(:,1),unique_spirals(:,2),unique_spirals_unit);
    xx = 1:1140; yy = 1:1320;
    [xxq,yyq] = meshgrid(xx,yy);
    vq1 = F(xxq,yyq);
    for i = 1:size(points,1)
        count_sample(i,1) = vq1(points(i,2),points(i,1));
    end
    count_sample(count_sample<0) = 0;
    %%
    subplot(2,7,7+count1);
    points1 = points-points(1,:);
    points_line = vecnorm(points1,2,2);
    plot(points_line,count_sample);
    xlim([0,800]);
    % ylim([0,4.5]);
    xticks([0 200 400 600 800]);
    xticklabels({'0','2','4','6','8'});
    ylabel('sprials/(mm^2*s)');    
    count1 = count1+1;
end
print(hs6cd, fullfile(save_folder,'Figs6cd_spirals_density_duration.pdf'),...
    '-dpdf', '-bestfit', '-painters');
