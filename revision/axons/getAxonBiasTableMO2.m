function getAxonBiasTableMO2(data_folder,save_folder)
%% get axon arbor bias for each neuron in sensory region with SVD
% first, define local axon regions, so that long range projections are excluded;
% then, do svd on axon terminal positions for each neuron.

% if soma is in VIS or RSP, no restriction;
% if soma is in auditory cortex, or secondary SS, restrict to within the
% area, because of cross-modality projections;
% if soma is in somatosensory cortex, restrict to the large SSp region, 
% so that long-range motor projections are excluded.
%% load singel cell morphology dataset
load(fullfile(data_folder,'axons','all_cell_with_parents.mat'));
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
av = readNPY(fullfile(data_folder,'tables',...
    'annotation_volume_10um_by_index.npy'));
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%% define axon arbor search regions by its soma location
MO_label = {'MOp','MOs','ACA','PL','ILA','ORB'};
MO_indx = getRegionIndex(MO_label,st);                           % get all index values in the cortical sensory area
for k = 1:6
    axon_indx{k} = MO_indx;
end
all_label = {'MOp','MOs','ACA','PL','ILA','ORB'};
%% get axon terminals for each cell, grouped by sensory regions
hemi = 1;                                                                  % only look at left hemisphere
for k = 1:numel(all_label)
    iregion = all_label(k);
    st_region_indx = getRegionIndex(iregion,st);                           % get brain region index
    [cell_id] = getCellPos(allCoords,av,st_region_indx);                   % get soma position in 3d atlas space
    if not(isempty(cell_id))
        [soma_all2{k},axon_terminal_all2{k}] = getAxonTerminal(...
            cell_id,allCoords,av,axon_indx{k}, hemi);                          % get axon terminals in specifed region boundries
        labels{k} = repmat(string(iregion),numel(soma_all2{k}),1);             % give an area label for each soma
    else
        soma_all2{k} = [];
        axon_terminal_all2{k} = [];
        lables{k} = [];
    end       
end
%% concatenate all cells
soma_all3 = cat(2, soma_all2{:});
axon_terminal_all3 = cat(2,axon_terminal_all2{:});
labels_all = cat(1,labels{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
labels_all = labels_all(a);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
%% main algorithm for axon terminal svd
[soma_center,axon_vector,angle1,polarity] = getAxonSVD(...
    axon_terminal_all3,soma_all3);   
%% calculate bias angle for each cell, reference to SSp center
center = [377,428]; % MOp
orthog_vector = soma_center-center;
for i =1:size(axon_vector,1)
    u = [axon_vector(i,:) 0];
    v = [orthog_vector(i,:) 0];
    ang_diff1(i,1) = atan2(norm(cross(u,v)),dot(u,v));
end    
ang_diff2 = rad2deg(ang_diff1);
%% calculate soma polar angle reference to SSp center
r = 1;
soma_angle = atan2(orthog_vector(:,2),orthog_vector(:,1));
z = r*exp(1i*soma_angle);
%% sort neuron by soma position
T = table(soma_center,soma_angle,z,axon_vector,...
    angle1,polarity,ang_diff2,labels_all);
T = renamevars(T,["z","axon_vector","angle1","polarity","ang_diff2"], ...
                 ["soma_polar_angle","axon_bias","bias_angle",...
                 "pc_ratio","center_bias_angle"]);
%%
writetable(T,fullfile(save_folder,'Axon_bias_all_cells_MO2.csv'));