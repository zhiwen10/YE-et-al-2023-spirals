function [soma_all1,axon_points_all1,dendrite_all1] = get_axon_points_by_neuron(cell_id,allCoords,av,axon_st_cortex_indx2, hemi)
type = 2; % only get axon
type1 = 1;
type2 = 4; % only get dendrite
% hemi =2; % both hemisphere
axon_points_all1 = {};
soma_all1 = {};
for i = 1:numel(cell_id)
    clear noi axonCoords axonCoords_filt
    axonCoords_filt = [];
    noi = allCoords{cell_id(i)}; % neuron of interest: brain ID 210254, neuron ID *20600
    axonCoords1 = noi(noi(:,1)== type,:);   
    % filter axon terminal in cortex
    [axonCoords_filt] = filter_axon_position(axonCoords1,av,axon_st_cortex_indx2,hemi);
    axon_points_all1{i} = axonCoords_filt;
    %%
    soma = noi(noi(:,1)== type1,2:4)/10;  
    soma_all1{i} = soma;
    %%
    dendrite = noi(noi(:,1)== type2,2:4)/10;  
    dendrite_all1{i} = dendrite;
end
end