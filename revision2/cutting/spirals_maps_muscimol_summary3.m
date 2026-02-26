githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'YE-et-al-2023-spirals')));            % paper repository
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data'; 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
spiral_folder = 'E:\task2\spirals\spirals_grouping';
data_folder = 'E:\task2';
T1 = readtable('cutting_session_list.xlsx');
%%
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
[sessions,ia,ic] = unique(T2.session_id);
T3 = T2(ia,:);
%%
rois = cell(numel(sessions),2);
for m = 1:numel(sessions)
    %%
    mn = T3.MouseID{m};
    session = sessions(m);    
    Ti = T1(ismember(T1.MouseID, mn) & T1.session_id == session & ismember(T1.label, 'cutting'),:);

    tda = Ti.date;
    en = Ti.folder;
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];

    load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
    session_root = fullfile(data_folder,'task_svd',fname);        
    load(fullfile(session_root,'meanImage.mat'));
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
    % [centers(:,:,m),radius(:,m)] = getCraniotomyROI(mimgt);
    imagesc(mimgt);
    roi1 = drawpolygon;
    rois{m,1} = roi1.Position;
    roi2 = drawpolygon;
    rois{m,2} = roi2.Position;
end
%%
save('craniotomy_registered2.mat','rois');
%%
label_all = {'spontaneous','cutting'};
spirals = cell(2,1);
frames = [];
for k = 1:2
    spirals_all = [];
    frameN_all = 0;
    for m = 1:numel(sessions)
        %%
        mn = T3.MouseID{m};
        session = sessions(m);    
        %%
        clear archiveCell filteredSpirals spiralsT lia
        label = label_all{k};
        Ti = T1(ismember(T1.MouseID, mn) & T1.session_id == session & ismember(T1.label, label),:);

        tda = Ti.date;
        en = Ti.folder;
        tdb = datestr(tda,'yyyymmdd');
        fname = [mn '_' tdb '_' num2str(en)];
        load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));
        load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
        session_root = fullfile(data_folder,'task_svd',fname);

        load(fullfile(session_root,'svdTemporalComponents_corr_timestamps.mat')); 
        load(fullfile(session_root,'meanImage.mat'));
        mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
        frameN = numel(t);

        spiral_duration = cellfun(@(x) size(x,1), archiveCell);
        indx2 = (spiral_duration>=2);
        groupedCells = archiveCell(indx2);
        filteredSpirals = cell2mat(groupedCells);
        [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,filteredSpirals(:,1),filteredSpirals(:,2));
        filteredSpirals(:,1:2) = round(spiralsT); 
        [lia,locb] = ismember(filteredSpirals(:,1:2),brain_index,'rows');
        filteredSpirals = filteredSpirals(lia,:);
        spirals_all = [spirals_all;filteredSpirals];
        frameN_all = frameN_all+frameN;            
    end
    spirals{k,1} = spirals_all;
    frames(k,1) = frameN_all;
end
%%
load('craniotomy_registered.mat');
th2 = 1:5:360; 
lineColor = 'w'; lineColor1 = 'w';
hemi = [];
scale3 = 5;
scale2 = 1;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
color1 = cbrewer2('qual','Set1',5);
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 800]);
ax0 = subplot(1,3,1);
m = 1;
mn = T3.MouseID{m};
session = sessions(m);    
label = label_all{k};
Ti = T1(ismember(T1.MouseID, mn) & T1.session_id == session & ismember(T1.label, 'cutting'),:);
tda = Ti.date;
en = Ti.folder;
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
session_root = fullfile(data_folder,'task_svd',fname);        
load(fullfile(session_root,'meanImage.mat'));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));

obj1 = imagesc(mimgt);
hold on;
scale2=1;
overlayOutlines(coords,scale2);
set(gca,'Ydir','reverse')
colormap(ax0,gray);
set(obj1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
axis off; axis image;
hold on;
for i = 1:numel(sessions)
    for j = 1:2
        clear pos
        pos = rois{i,j};
        pos(end+1,:) = pos(1,:);
        plot(pos(:,1),pos(:,2),'color',color1(i,:));
        hold on;
    end
end
% for id = 1:5
%     for j = 1:2
%         px1 = centers(j,1,id);
%         py1 = centers(j,2,id);
%         r = radius(j,id);
%         cx2 = r*cosd(th2)+px1;
%         cy2 = r*sind(th2)+py1;
%         hold on;
%         plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1(id,:),'LineWidth',1);          % draw the circle at max radius
%     end
% end
for k = 1:2
    clear sprials_temp unique_spirals unique_spirals_unit
    ax1 = subplot(1,3,k+1);
    sprials_temp = spirals{k,1};
    frame_all = frames(k,1);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(sprials_temp,hist_bin);
    unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
    [ax1,cb1]= plotDesnity(ax1,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
    caxis([0,2]);
end
%%
print(h1, 'SSp_cutting_summary4', '-dpdf', '-bestfit', '-painters');

