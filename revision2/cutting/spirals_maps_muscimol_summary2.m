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
for m = 1:numel(sessions)
    %%
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
    [centers(:,:,m),radius(:,m)] = getCraniotomyROI(mimgt);
end
save('craniotomy_registered.mat','centers','radius');
%%
load('craniotomy_registered.mat');
th2 = 1:5:360; 
lineColor = 'w'; lineColor1 = 'w';
hemi = [];
scale3 = 5;
scale2 = 1;
label_all = {'spontaneous','cutting'};
h1 = figure('Renderer', 'painters', 'Position', [100 100 900 800]);
for m = 1:numel(sessions)
    %%
    mn = T3.MouseID{m};
    session = sessions(m);    
    %%
    for k = 1:2
        %%
        clear archiveCell filteredSpirals spiralsT lia
        %%
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
        frame_all = numel(t);
        frameLim = 0;
        ax1 = subplot(5,3,(m-1)*3+k+1);
        [unique_spirals,unique_spirals_unit] = getSpiralDensity3(archiveCell,tform,brain_index,frameLim,frame_all);
        [ax1,cb1]= plotDesnity(ax1,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
        caxis([0,2]);
        
        if k == 2
            ax0 = subplot(5,3,(m-1)*3+1);
            obj1 = imagesc(mimgt);
            hold on;
            overlayOutlines(coords,scale2);
            set(gca,'Ydir','reverse')
            colormap(ax0,gray);
            set(obj1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
            hold on;
            plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
            axis off; axis image;
            cb2 = colorbar;
            cb2.Position(1) = cb2.Position(1)+0.05;
            cb2.Position(2) = cb2.Position(2)+0.1;
            cb2.Position(4) = cb2.Position(4)/2;
            hold on;
            for j = 1:2
                px1 = centers(j,1,m);
                py1 = centers(j,2,m);
                r = radius(j,m);
                cx2 = r*cosd(th2)+px1;
                cy2 = r*sind(th2)+py1;
                hold on;
                plot([cx2 cx2(1)],[cy2 cy2(1)],'w','LineWidth',1);          % draw the circle at max radius
            end
        end
            
    end
end
%%
print(h1, 'SSp_cutting_summary2', '-dpdf', '-bestfit', '-painters');

