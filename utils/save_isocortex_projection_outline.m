data_folder = 'C:\Users\Steinmetz lab\Documents\git\YE-et-al-2023-spirals-test\data';
[projectedAtlas1,projectedTemplate1,coords] = get_isocortex_horizontal_projection(data_folder);
save(fullfile(data_folder,'isocortex_horizontal_projection_outline.mat'),'projectedAtlas1','projectedTemplate1','coords');
%%
function [projectedAtlas1,projectedTemplate1,coords] = get_isocortex_horizontal_projection(data_folder)
%% load 10um horizontal atlas and outline
load(fullfile(data_folder, 'ctxOutlines.mat'));
load(fullfile(data_folder,'projectedOutlineAtlas.mat'))
st = loadStructureTree(fullfile(data_folder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas1(projectedAtlas,projectedTemplate,st);
end