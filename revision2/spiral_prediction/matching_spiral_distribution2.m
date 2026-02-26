%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
load(fullfile(data_folder,'spirals','spirals_symmetry','roiSelection.mat'));
%%
MatchPointS_all = [];
for kk = 1:15
    %%
    clear spiralsT tf2 tf4 pointSL pointSR MatchPointSL MatchPointSR filteredSpirals2
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals2 = cell2mat(groupedCells); 
    %%
    filteredSpirals2(filteredSpirals2(:,4) == -1,4) = 0;
    spirals2 = filteredSpirals2;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,spirals2(:,1),spirals2(:,2));
    spirals2(:,1:2) = spiralsT;    
    tf2 = inROI(roiSL,spirals2(:,1),spirals2(:,2)); 
    pointSL = spirals2(tf2,:);
    tf4 = inROI(roiSR,spirals2(:,1),spirals2(:,2)); 
    pointSR = spirals2(tf4,:);
    %%
    [MatchPointS] = UniqueMatchPoints2(pointSL,pointSR);   
    %%
    MatchPointS_all = [MatchPointS_all;MatchPointS];
end
%%
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
load(fullfile(data_folder,'spirals','spirals_density',...
    'histogram_40pixels.mat'));
%%
hist_bin = 40;
[unique_spirals,scolor,low_color_bound,high_color_bound] = ...
    density_color_plot(MatchPointS_all,hist_bin);
%%
h1e = figure;
hist_bin = 40;
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);     % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_total*35;                 % spirals/(mm^2*s)
ax1 = subplot(1,1,1);
scatter(unique_spirals(:,1),unique_spirals(:,2),...
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