function getSpiralsDensityLine_newsession(data_folder,save_folder)
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
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
for kk = 1:4
    %%
    clear spiralsT filteredSpirals unique_spirals unique_spirals_unit indx2
    mn = fnames{kk};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T = T_session(T_session.label == "passive",:);
    %%
    session = 5;
    mn = T.MouseID{session};
    tda = T.date(session);
    en = T.folder(session);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    t = readNPY(fullfile(data_folder,'task','task_svd',fname,...
        'svdTemporalComponents_corr.timestamps.npy'));                     % read how many total frames
    nframe = numel(t);  
    %%
    load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
    load(fullfile(data_folder,'task','spirals_all','spirals_grouping',...
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
save(fullfile(save_folder,'spiralDensityLinePerSession_new.mat'),...
    'count_sample');