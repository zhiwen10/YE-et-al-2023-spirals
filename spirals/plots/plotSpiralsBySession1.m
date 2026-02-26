function hs5a1 = plotSpiralsBySession1(T,session_rows,data_folder,save_folder)
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
hs5a1 = figure('Renderer', 'painters', 'Position', [100 100 1100 450]);
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
    %% plot receptive field sign map
    ax0(count1) = subplot(1,6,count1);                                     % plot registered sign map with brain outline
    im = imagesc(ax0(count1),signMaptransformed);                          % receptive field sign map
    set(ax0(count1), 'Color', 'w');
    set(ax0(count1), 'XTickLabel', '', 'YTickLabel', '');
    axis image; axis off;
    colormap(ax0(count1),colormap_RedWhiteBlue)
    set(im, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
    hold on;
    overlayOutlines(coords,scale,'k');
    cb00 = colorbar(ax0(count1));
    cb00.Ticks = [-1,0,1];
    cb00.TickLabels = {'-1','0','1'};
    cb_size00 =  get(cb00,'Position');                                     % gets the positon and size of the color bar
    set(cb00,'Position',...
        [cb_size00(1)+0.04 cb_size00(2)+0.03 ...
        cb_size00(3) cb_size00(4)/2]);                                     % To adjust plot position
    fname2 = [mn '_' tdb '_' num2str(en)];
    title(fname2,'Interpreter','None');     
    %%
    count1 = count1+1;
end
fig_name = ['FigS5a1_rfmap_by_session_' num2str(session_rows(1)) '_' ,...
    num2str(session_rows(end))];
print(hs5a1,fullfile(save_folder,fig_name),...
    '-dpdf', '-bestfit', '-painters');                                     % save as pdf file