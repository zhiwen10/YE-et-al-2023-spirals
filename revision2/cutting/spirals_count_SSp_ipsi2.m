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
spath = string(st.structure_id_path);
%%
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%% right SSp and MO index
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
hemi = [];
scale = 1;
[indexSSp,UselectedSSp] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
hemi = [];
[indexMO,UselectedMO] = select_area(...
    frotalArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MO_L = BW_empty; BW_MO_L(indexMO) = 1;
BW_SSp_L = BW_empty; BW_SSp_L(indexSSp) = 1;
BW_MO_L(:,570:end) = 0; BW_SSp_L(:,570:end) = 0;
[row1,col1] = find(BW_MO_L); MO_L_index = [col1,row1];
[row2,col2] = find(BW_SSp_L); SSp_L_index = [col2,row2];

BW_MO_R = BW_empty; BW_MO_R(indexMO) = 1;
BW_SSp_R = BW_empty; BW_SSp_R(indexSSp) = 1;
BW_MO_R(:,1:570) = 0; BW_SSp_R(:,1:570) = 0;
[row1,col1] = find(BW_MO_R); MO_R_index = [col1,row1];
[row2,col2] = find(BW_SSp_R); SSp_R_index = [col2,row2];

area_index = {MO_L_index,MO_R_index,SSp_L_index,SSp_R_index};
%%
spiral_folder = 'E:\task2\spirals\spirals_grouping';
data_folder = 'E:\task2';
T1 = readtable('cutting_session_list.xlsx');
%%
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
%%
spiral_density = [];
for m = 1:size(T2,1)
    %%
    clear archiveCell filteredSpirals spiralsT lia
    %%
    mn = T2.MouseID{m};
    tda = T2.date(m);
    en = T2.folder(m);
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));
    load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
    session_root = fullfile(data_folder,'task_svd',fname);
    load(fullfile(session_root,'svdTemporalComponents_corr_timestamps.mat')); %t
    %%
    frameLim = 0;
    frame_all = numel(t);
    %%
    spiral_duration = cellfun(@(x) size(x,1), archiveCell);
    indx2 = (spiral_duration>=2);
    groupedCells = archiveCell(indx2);
    filteredSpirals = cell2mat(groupedCells);
    filteredSpirals = filteredSpirals(filteredSpirals(:,5)>frameLim,:);
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,filteredSpirals(:,1),filteredSpirals(:,2));
    filteredSpirals(:,1:2) = round(spiralsT); 
    for kk = 1:4
        clear lia locb spirals
        [lia,locb] = ismember(filteredSpirals(:,1:2),area_index{kk},'rows');
        spirals = filteredSpirals(lia,:);
        spiral_count = size(spirals,1);
        spiral_density(m,kk) = spiral_count./frame_all*35;             
    end
end
%%
Ta = array2table(spiral_density,'VariableNames',{'MOL','MOR','SSpL','SSpR'});
Ta1 = [T2,Ta];
%%
sessions = unique(Ta1.session_id);
MO_list = []; SSp_list = [];
for i = 1:numel(sessions)
    clear Ta2 MO_list1 SSp_list1 
    Ta2 = Ta1(Ta1.session_id == sessions(i),:);
    MO_list1 = [Ta2.MOL';Ta2.MOR'];
    SSp_list1 = [Ta2.SSpL';Ta2.SSpR'];
    MO_list = [MO_list;MO_list1];
    SSp_list = [SSp_list;SSp_list1];
end
%%
clear MO_list2 SSp_list2
for i = 1:numel(sessions)
    MO_list2(i,:) = MO_list(2*i-1,:)+ MO_list(2*i,:);
    SSp_list2(i,:) = SSp_list(2*i-1,:)+ SSp_list(2*i,:);
end
MO_list = MO_list2;
SSp_list = SSp_list2;
%%
MO_list = MO_list(:,[1,3]);
SSp_list = SSp_list(:,[1,3]);
%%
spiral_density_control2 = [MO_list(:,1),SSp_list(:,1)];
spiral_density_craniotomy2 = [MO_list(:,2),SSp_list(:,2)];
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 700 400]);
for i = 1:2
    scatter(ones(size(spiral_density_control2))*(2*i-1),spiral_density_control2(:,i),8,'k');
    hold on;
end
for i = 1:2
    scatter(ones(size(spiral_density_control2))*2*i,spiral_density_craniotomy2(:,i),8,'r');
    data_temp1 = spiral_density_control2(:,i);
    data_temp2 = spiral_density_craniotomy2(:,i);
    plot([ones(size(data_temp1))*(2*i-1),ones(size(data_temp2))*(2*i)]',...
        [data_temp1,data_temp2]','color','k');
end
xticks([1.5:2:3.5]);
xticklabels({'MO','SSp'});
ylim([0,20]);
%%
[h1,p1] = ttest(MO_list(:,1),MO_list(:,2));
%%
print(h1, 'Cutting_SSp_summary_scatter2', '-dpdf', '-bestfit', '-painters');    