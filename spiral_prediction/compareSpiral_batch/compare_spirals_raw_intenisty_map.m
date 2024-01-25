% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%% load 10um horizontal atlas and outline
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF-Zhiwen\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2); st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas,idFilt);
projectedAtlas1(~Lia) = 0; projectedTemplate1(~Lia) = 0;
areaPath{1} = '/997/8/567/688/695/315/453/322/'; % SSp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
scale = 1;
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,'right',scale);
[indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,'left',scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_right = BW_empty; BW_right(indexright) =1;
BW_left = BW_empty; BW_left(indexleft) =1;
%%
[rowL,colL] = find(BW_left);
brain_index_left = [colL,rowL];
[rowR,colR] = find(BW_right);
brain_index_right = [colR,rowR];
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T1 = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
raw_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_raw_fftn';
predict_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_predict';
reg_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\registration';
%%
spiral_left_match_all = {};
spiral_right_match_all = {};
ratio_all = [];
ratio =[];
for kk = 1:size(T1,1)
    clear spiralsT spirals_filt spirals_prediction_left spirals_prediction_right spiralsT1 sprials_left spirals_right
    ops = get_session_info(T1,kk);
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    serverRoot = expPath(ops.mn, ops.td, ops.en);   
    t = readNPY(fullfile(serverRoot, 'blue','svdTemporalComponents.timestamps.npy'));
    load(fullfile(raw_folder,[fname '_spirals_group_fftn.mat']));
    load(fullfile(reg_folder,[fname '_tform.mat']));
    spiral_length = cellfun(@(x) size(x,1), archiveCell);
    spiral_sequence = archiveCell(spiral_length>=2);
    spirals_filt= cell2mat(spiral_sequence);
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,spirals_filt(:,1),spirals_filt(:,2));
    spirals_filt(:,1:2) = round(spiralsT); 
    %%
    [liaL,locbL] = ismember(spirals_filt(:,1:2),brain_index_left,'rows');
    [liaR,locbR] = ismember(spirals_filt(:,1:2),brain_index_right,'rows');
    spirals_left = spirals_filt(liaL,:);
    spirals_right = spirals_filt(liaR,:);
    %%
    frame_all = numel(t);
    clear pwAll
    %%
    hist_bin = 40;
    [unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spirals_filt,hist_bin);
    pixSize = 0.01; % mm/pix
    pixArea = pixSize^2;
    h = figure;
    unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
    unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
    scatter(unique_spirals(:,1),unique_spirals(:,2),3,unique_spirals_unit,'filled');
    hold on;
    scale2=1;
    overlayOutlines(coords,scale2);
    set(gca,'Ydir','reverse')
    colormap(hot);
    axis image;
    colorbar;
    title(fname);
end