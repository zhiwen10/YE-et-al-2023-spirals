function hs5a2 = plotSpiralsBySession2(T,session_rows,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));
BW = logical(projectedAtlas1);                                             % atlas brain boundry binary mask 
[row,col] = find(BW);
brain_index = [col,row];                                                   % position idenx within the brain boundry
%% set params
hist_bin = 40;                                                             % bin size (pixels) for spiral density estimation
pixSize = 0.01;                                                            % pixel resolution mm/pix
pixArea = pixSize^2;                                                       % 2d-pixel resolution mm^2/pix^2
scale  = 1;                                                                % no need to scale for atlas outline here
%%
hs5a2 = figure('Renderer', 'painters', 'Position', [100 100 1100 450]);
count1 = 1;
for kk = session_rows                                                      % plot 6 session in one figure, to avoid crowding
    %% session info
    clear spiralsT 
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load data
    fname = [mn '_' tdb '_' num2str(en)];
    t = readNPY(fullfile(data_folder,'spirals','svd',fname,...
        'svdTemporalComponents_corr.timestamps.npy'));                     % read how many total frames
    frame_all = numel(t);    
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));                               % load spiral centers (>40 pixels radius)
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));           % load atlas transformation matrix tform
    load(fullfile(data_folder,'spirals','rf_signmap',[fname '_signMap.mat']));       % load receptive field signMap
    mimgtransformed = imwarp(mimg,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));                   % transform mean wf image to atlas space
    signMaptransformed = imwarp(signMap,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));                   % transform signMap image to atlas space
    %% only use spiral sequences with at least 2 consecutive frames
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);                                          % sprial length > = 2 frames
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));                  % transform sprials to atlas space
    filteredSpirals(:,1:2) = round(spiralsT);    
    [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');      % find spirals within the brain boundry
    filteredSpirals = filteredSpirals(lia,:);    
    %% plot spirals density map
    ax1(count1) = subplot(1,6,count1);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(filteredSpirals,hist_bin);                      % hsitogram based density estimate with red color scale
    unique_spirals_unit = ...
        unique_spirals(:,3)/(hist_bin*hist_bin*pixArea);                   % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./frame_all(1)*35;            % spirals/(mm^2*s)
    scatter(ax1(count1),unique_spirals(:,1),unique_spirals(:,2),...
        1,unique_spirals_unit,'filled');                                   % plot scatter plot with color indicating spiral density 
    max_c = 1.0;
    caxis([0,max_c]);
    colormap(ax1(count1),hot);
    cb0 = colorbar(ax1(count1)); 
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    axis off; axis image;
    xlim(ax1(count1),[0,1140]); ylim(ax1(count1),[0,1320]);
    
    cb0.Ticks = [0,max_c];
    cb0.TickLabels = {'0',num2str(max_c)};
    cb_size0 =  get(cb0,'Position');                                       % gets the positon and size of the color bar
    set(cb0,'Position',...
        [cb_size0(1)+0.04 cb_size0(2)+0.03 ...
        cb_size0(3) cb_size0(4)/2]);                                       % To adjust plot position      
    %%
    count1 = count1+1;
end
fig_name = ['FigS5a2_spirals_by_session_' num2str(session_rows(1)) '_' ,...
    num2str(session_rows(end))];
print(hs5a2,fullfile(save_folder,fig_name),...
    '-dpdf', '-bestfit', '-painters');                                     % save as pdf file