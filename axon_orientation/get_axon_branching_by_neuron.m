function [soma_all,axon_all] = get_axon_branching_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx, hemi)
%%
type = 2; % only get axon
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
    [axonCoords_filt] = filter_axon_position(axonCoords,av,axon_st_cortex_indx, hemi);
    axon_all{i} = axonCoords_filt;
end
%% get soma
type1 = 1;
for i = 1:numel(cell_id)
    clear noi somaCoords
    somaCoords = [];
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    somaCoords = noi(noi(:,1)==type1,2:4)/10;
    soma_all{i} = somaCoords;
end