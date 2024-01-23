function [soma_all,axon_all] = get_all_axon_in_area_by_neuron(st,av,allCoords,region_label,axon_st_cortex_indx,cell_limit)
%%
% region_label = {"SSp-ll"};
st_region_indx = get_region_index_list(region_label,st); % get all index values that contain "SSp-ll" in st table
[cell_id] = get_cell_id(allCoords,av,st_region_indx); % get all cell index that have cell body in "SSp-ll"
%%
if nargin==5 % if cell_limit is not specified 
    cell_limit = cell_id;
else 
    if cell_limit>numel(cell_id)
        cell_limit = numel(cell_id);
    end
    cell_id_index = randperm(numel(cell_id),cell_limit);
    cell_id = cell_id(cell_id_index);
end
%%
type = 2; % only get axon
hemi = 1; % only get axons in left hemispehre
for i = 1:numel(cell_id)
    clear noi axonCoords axonCoords_filt
    axonCoords_filt = [];
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    parentsa = unique(noi(:,6)); %find parents
    x__n = noi(:,5);
    [x1,x2]=ismember(x__n,parentsa);
    noi1 = noi(not(x1),:);
    axonCoords = noi1(noi1(:,1)== type,:);   
    % filter axon terminal in cortex
    % axon_st_cortex_indx = get_region_index_list(region_label,st); % get all index values that contain "Isocortex" in st table
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx, hemi);
    axon_all{i} = axonCoords_filt;
end
% get soma
type1 = 1;
for i = 1:numel(cell_id)
    clear noi somaCoords
    somaCoords = [];
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    somaCoords = noi(noi(:,1)==type1,2:4)/10;
    soma_all{i} = somaCoords;
end