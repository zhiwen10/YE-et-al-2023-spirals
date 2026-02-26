function getSpiralsDensityLine(T,data_folder,save_folder)
%%
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
pixSize = 0.01;                                                            % mm/pix after registration
pixArea = pixSize^2;
%% draw a line 
x = [105,520]; y = [780,225];
xq = 105:520;
vq = interp1(x,y,xq);
points = [xq',round(vq)'];
xx = 1:1140; yy = 1:1320;
[xxq,yyq] = meshgrid(xx,yy);
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%   
for kk = 1:15
    %%
    clear spiralsT filteredSpirals unique_spirals unique_spirals_unit indx2
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    t = readNPY(fullfile(data_folder,'spirals','svd',fname,...
        'svdTemporalComponents_corr.timestamps.npy'));                     % read how many total frames
    nframe = numel(t);  
    %%
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);   
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,filteredSpirals(:,1),filteredSpirals(:,2));    
    filteredSpirals(:,1:2) = round(spiralsT); 
    hist_bin = 40;
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(filteredSpirals,hist_bin);
    %%
    cell_count = size(filteredSpirals,1);
    unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./nframe*35;                  % spirals/(mm^2*s)
    % interp histgram counts 
    F = scatteredInterpolant(unique_spirals(:,1),unique_spirals(:,2),...
        unique_spirals_unit);
    vq1 = F(xxq,yyq);
    for i = 1:size(points,1)
        count_sample(i,kk) = vq1(points(i,2),points(i,1));
    end
end
save(fullfile(save_folder,'spiralDensityLinePerSession.mat'),...
    'count_sample');