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
% areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
% areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
% areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
% areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
% areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
% areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
% areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
% areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
% areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
hemi = [];
scale = 1;
[indexSSp,UselectedSSp] = select_area(...
    sensoryArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);

% MOp index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/985/"; % MOp
MOPArea = strcat(areaPath(:));
hemi = [];
[indexMOP,UselectedMOP] = select_area(...
    MOPArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);

% MOs index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/993/"; % MOs
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
MOSArea = strcat(areaPath(:));
hemi = [];
[indexMOS,UselectedMOS] = select_area(...
    MOSArea,spath,st,coords,Utransformed,...
    projectedAtlas1,projectedTemplate1,hemi,scale);
%%
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_MOP_L = BW_empty; BW_MOP_L(indexMOP) = 1;
BW_MOS_L = BW_empty; BW_MOS_L(indexMOS) = 1;
BW_SSp_L = BW_empty; BW_SSp_L(indexSSp) = 1;
BW_MOP_L(:,570:end) = 0; BW_MOS_L(:,570:end) = 0; BW_SSp_L(:,570:end) = 0;
[row1,col1] = find(BW_MOP_L); MOP_L_index = [col1,row1];
[row2,col2] = find(BW_SSp_L); SSp_L_index = [col2,row2];
[row3,col3] = find(BW_MOS_L); MOS_L_index = [col3,row3];

BW_MOP_R = BW_empty; BW_MOP_R(indexMOP) = 1;
BW_MOS_R = BW_empty; BW_MOS_R(indexMOS) = 1;
BW_SSp_R = BW_empty; BW_SSp_R(indexSSp) = 1;
BW_MOP_R(:,1:570) = 0; BW_MOS_R(:,1:570) = 0; BW_SSp_R(:,1:570) = 0;
[row1,col1] = find(BW_MOP_R);  MOP_R_index = [col1,row1];
[row2,col2] = find(BW_SSp_R); SSp_R_index = [col2,row2];
[row3,col3] = find(BW_MOS_R); MOS_R_index = [col3,row3];

area_index = {MOS_L_index,MOS_R_index,MOP_L_index,MOP_R_index,SSp_L_index,SSp_R_index};
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
    for kk = 1:6
        clear lia locb spirals
        [lia,locb] = ismember(filteredSpirals(:,1:2),area_index{kk},'rows');
        spirals = filteredSpirals(lia,:);
        spiral_count = size(spirals,1);
        spiral_density(m,kk) = spiral_count./frame_all*35;             
    end
end
%%
Ta = array2table(spiral_density,'VariableNames',{'MOSL','MOSR','MOPL','MOPR','SSpL','SSpR'});
Ta1 = [T2,Ta];
Ta1.MOS = Ta1.MOSL+Ta1.MOSR;
Ta1.MOP = Ta1.MOPL+Ta1.MOPR;
Ta1.SSp = Ta1.SSpL+Ta1.SSpR;
spiral_density_control2 = Ta1{1:3:end,[15,16,17]};
spiral_density_cutting2 = Ta1{3:3:end,[15,16,17]};
%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 700 400]);
for i = 1:3
    scatter(ones(size(spiral_density_control2))*(2*i-1),spiral_density_control2(:,i),8,'k');
    hold on;
    scatter(ones(size(spiral_density_control2))*2*i,spiral_density_cutting2(:,i),8,'r');
end
for i = 1:3    
    data_temp1 = spiral_density_control2(:,i);
    data_temp2 = spiral_density_cutting2(:,i);
    plot([ones(size(data_temp1))*(2*i-1),ones(size(data_temp2))*(2*i)]',...
        [data_temp1,data_temp2]','color','k');
end
xticks([1.5:2:7.5]);
xticklabels({'MOs','MOp','SSp'});
ylim([0,20]);
%%
[h1,p1] = ttest(spiral_density_control2(:,1),spiral_density_cutting2(:,1));
[h2,p2] = ttest(spiral_density_control2(:,2),spiral_density_cutting2(:,2));
%%
print(h1, 'Cutting_SSp_summary_scatter3', '-dpdf', '-bestfit', '-painters');    