function [soma_all,axon_all] = getAxonTerminal(cell_id,allCoords,av,filter_region, hemi)
%%
% morphology dataset rows are organized as below:
% pointType,x,y,z,unique_point_id, point_parent_id 
% 3,3159,3841,3337,210,209;
% 3,3162,3845,3337,211,210;
% 3,3127,3676,3330,212,81;
% 3,3123,3676,3328,213,212;
% 3,3119,3677,3327,214,213;
% 3,3115,3677,3326,215,214

% pointType indicates if the point is a soma [1], an axon [2], or an
% dendrite [3]

% x,y,z positions are registered to allen atlas with 1um resolution
% each point has a unique id [210] in column 5
% each point has a parent id [209] in column 6, which indicates what its
% upstream point is

% if a unique point is a terminal, then no other points use it as parent
% therefore, its unique index won't show up in the parent column
%% get axon terminal 3d coordinates
type = 2; 
for i = 1:numel(cell_id)
    clear noi axonCoords axonCoords_filt
    axonCoords_filt = [];
    noi = allCoords{cell_id(i)}; 
    parentsa = unique(noi(:,6));                                           % find all point id that had a downstream point [column6]
    x_n = noi(:,5);
    [x1,x2]=ismember(x_n,parentsa);                                       % find which unique ids [column5] were parents
    noi1 = noi(not(x1),:);                                                 % find unique ids that were not parents, and are therefore terminals
    axonCoords = noi1(noi1(:,1)== type,:);                                 % terminals has to be type 2 in [column1]
    [axonCoords_filt] = filter_axon_position(...
        axonCoords,av,filter_region, hemi);                                % filter axon terminal in specified region and hemi, in 10um resolution
    axon_all{i} = axonCoords_filt;
end
%% get soma 3d coordinates
type1 = 1;
for i = 1:numel(cell_id)
    clear noi somaCoords
    somaCoords = [];
    noi = allCoords{cell_id(i)}; 
    somaCoords = noi(noi(:,1)==type1,2:4)/10;                              % find soma [type1] 3d coordinates, in 10um resolution
    soma_all{i} = somaCoords;
end