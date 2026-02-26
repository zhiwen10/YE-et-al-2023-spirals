function hs7a = plotExampleSpiralTrajectory(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[areaPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
[Lia,Locb] = ismember(projectedAtlas1,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
sensoryArea = strcat(areaPath(:));

spath3 = startsWith(spath,sensoryArea);
idFilt1 = st.index(spath3);
[Lia1,Locb1] = ismember(projectedAtlas1,idFilt1);
projectedAtlas1(~Lia1) = 0; 
projectedTemplate1(~Lia1) = 0;
hemi = 'right';
if strcmp(hemi,'right')
    projectedAtlas1(:,1:size(projectedAtlas1,2)/2) = 0; 
    projectedTemplate1(:,1:size(projectedTemplate1,2)/2) = 0;
elseif strcmp(hemi,'left')
    projectedAtlas1(:,size(projectedAtlas1,2)/2+1:end) = 0; 
    projectedTemplate1(:,size(projectedAtlas1,2)/2+1:end) = 0; 
end
BW_SSp = logical(projectedAtlas1);
[row1,col1] = find(BW_SSp);
brain_index1 = [col1,row1];
%% mask for whole brain 
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
kk = 7;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];

load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
load(fullfile(data_folder,'spirals','spirals_grouping',...
    [fname '_spirals_group_fftn.mat']));
%%
clear indx2 groupedCells
indx2 = cellfun(@(x) size(x,1), archiveCell);
indx3 = (indx2>=7);
groupedCells = archiveCell(indx3);
%%
colorn = max(indx2+5);
color_cw = flipud(cbrewer2('seq','Blues',colorn));
color_ccw = flipud(cbrewer2('seq','OrRd',colorn));
scale = 1;
%%
hs7a = figure;
overlayOutlines(coords,scale,[0.8,0.8,0.8]);
set(gca, 'YDir','reverse');
for i = 1:numel(groupedCells)
    clear pwAllEpoch1 u v
    pwAllEpoch1 = groupedCells{i};
    [u,v] = transformPointsForward(tform,pwAllEpoch1(:,1),pwAllEpoch1(:,2));
    hold on;
    if pwAllEpoch1(1,4) ==0
        cc = color_cw;
    else
        cc = color_ccw;
    end
    for k = 1:numel(u)-1
        plot([u(k),u(k+1)],[v(k),v(k+1)],'color',cc(k+5,:),'lineWidth',1.5);   
        hold on;
    end
    scatter(u(1),v(1),18,'MarkerFaceColor',cc(6,:),'MarkerEdgeColor','None');
    axis image; axis off;
end
%%
print(hs7a, fullfile(save_folder,'Figs7a_spirals_trajectory'),...
    '-dpdf', '-bestfit', '-painters');