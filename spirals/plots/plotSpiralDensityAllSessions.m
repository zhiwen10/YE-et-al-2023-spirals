function h1e = plotSpiralDensityAllSessions(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
load(fullfile(data_folder,'spirals','spirals_density',...
    'histogram_40pixels.mat'));
%%
hist_bin = 40;
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_total*35;                 % spirals/(mm^2*s)
unique_spirals(:,3) = unique_spirals_unit;
% save('spiral_density.mat','unique_spirals');
%%
h1e = figure;
ax1 = subplot(1,1,1);
scatter(unique_spirals(:,1),unique_spirals(:,2),...
    3,unique_spirals(:,3),'filled');
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
%%
print(h1e, fullfile(save_folder,'Fig1e_spiral_histogram'), ...
    '-dpdf', '-bestfit', '-painters');
%%
% Map data values to hot colormap
cmax = max(unique_spirals_unit(:));
n_colors = 256;
cmap = hot(n_colors);

% Normalize data to [0, 1] based on caxis range [0, cmax]
vals = unique_spirals(:,3);
vals_norm = vals / cmax;
vals_norm = max(0, min(1, vals_norm));  % clamp to [0, 1]

% Map to colormap indices
color_idx = round(vals_norm * (n_colors - 1)) + 1;
vertex_colors = cmap(color_idx, :);  % Nx3 RGB

% Create 1320x1140x3 image (rows=Y, cols=X)
density_colorrgb = zeros(1320, 1140, 3);

% Fill in pixel values from scatter coordinates
x = round(unique_spirals(:,1));
y = round(unique_spirals(:,2));

% Clamp to valid range
x = max(1, min(1140, x));
y = max(1, min(1320, y));

for i = 1:length(x)
    density_colorrgb(y(i), x(i), :) = vertex_colors(i, :);
end

% Save
% save( 'density_colorrgb.mat', 'density_colorrgb');