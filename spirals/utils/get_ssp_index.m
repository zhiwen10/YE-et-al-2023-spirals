
function brain_index = get_ssp_index(data_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));

% only select cortex in the atlas
st = loadStructureTree(fullfile(data_folder,'tables',...
    'structure_tree_safe_2017.csv'));                                      % a table of what all the labels mean
spath = string(st.structure_id_path);
spath2 = startsWith(spath,'/997/8/567/');
idFilt = st.index(spath2);
st1  = st(spath2,:);
[Lia,Locb] = ismember(projectedAtlas1,idFilt);
projectedAtlas1(~Lia) = 0; 
projectedTemplate1(~Lia) = 0;
areaPath{1} = '/997/8/567/688/695/315/453/322/'; % SSp
sensoryArea = strcat(areaPath(:));
Utransformed = zeros(size(projectedAtlas1));
hemi = 'right';
scale = 1;
[indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
BW_right = BW_empty; 
BW_right(indexright) =1;
[row,col] = find(BW_right);
brain_index = [col,row];