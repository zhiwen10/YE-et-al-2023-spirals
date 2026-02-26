function hs7b = plotSpiralDirectionRatio(T,data_folder,save_folder)
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
%% mask for whole brain 
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];

projectedAtlas_right = projectedAtlas1;
projectedAtlas_right(:,1:size(projectedAtlas_right,2)/2) = 0; 
BW_SSp_right = logical(projectedAtlas_right);
[row_right,col_right] = find(BW_SSp_right);
brain_index_right = [col_right,row_right];

projectedAtlas_left = projectedAtlas1;
projectedAtlas_left(:,size(projectedAtlas_left,2)/2+1:end) = 0; 
BW_SSp_left = logical(projectedAtlas_left);
[row_left,col_left] = find(BW_SSp_left);
brain_index_left = [col_left,row_left];
%%
for kk = 1:15
    clear spiralsT spiral2 spiral3 filteredSpirals2 spiralsT
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals2 = cell2mat(groupedCells);
    %%
    filteredSpirals2(filteredSpirals2(:,4) == -1,4) = 0;
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,filteredSpirals2(:,1),filteredSpirals2(:,2));
    filteredSpirals2(:,1:2) = spiralsT;    
    filteredSpirals2(:,1:2) = round(filteredSpirals2(:,1:2));
    %% mask for whole brain
    [lia,locb] = ismember(filteredSpirals2(:,1:2),brain_index,'rows');
    filteredSpirals2 = filteredSpirals2(lia,:);
    %%
    projectedAtlas_right = projectedAtlas1;
    projectedAtlas_right(:,1:size(projectedAtlas_right,2)/2) = 0; 
    BW_SSp_right = logical(projectedAtlas_right);
    [row_right,col_right] = find(BW_SSp_right);
    brain_index_right = [col_right,row_right];
    %%
    [lia_right,locb_right] = ismember(filteredSpirals2(:,1:2),brain_index_right,'rows');
    spirals_right = filteredSpirals2(lia_right,:);
    [lia_left,locb_left] = ismember(filteredSpirals2(:,1:2),brain_index_left,'rows');
    spirals_left = filteredSpirals2(lia_left,:);
    direction_ratio(kk,1) = sum(spirals_left(:,4))/size(spirals_left,1); % 1 is counterclockwise
    direction_ratio(kk,2) = sum(spirals_right(:,4))/size(spirals_right,1);
end
%% t-test
[ha,p] = ttest(direction_ratio(:,1),direction_ratio(:,2));
%%
mean_ratio = mean(direction_ratio);
std_ratio = std(direction_ratio);
hs7b = figure('Renderer', 'painters', 'Position', [100 100 200 250])
bar([0,1],mean_ratio,'FaceColor',[0.8,0.8 0.8],'EdgeColor',[0 0 0]);
hold on;
scatter(0,direction_ratio(:,1),36,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');
hold on;
scatter(1,direction_ratio(:,2),36,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k');
hold on;
for k = 1:size(direction_ratio,1)
    plot([0,1],[direction_ratio(k,1),direction_ratio(k,2)],'k');
    hold on;
end
errorbar(0,mean_ratio(:,1),std_ratio(:,1),'lineWidth',1,"CapSize",15,'color','k');
hold on;
errorbar(1,mean_ratio(:,2),std_ratio(:,2),'lineWidth',1,"CapSize",15,'color','k');
ylim([0,0.6]);
set(gca,'xticklabel',{'left','right'})
ylabel('ccw spiral ratio');
print(hs7b, fullfile(save_folder,'FigS7b_spirals_direction_all.pdf'), ...
    '-dpdf', '-bestfit', '-painters');