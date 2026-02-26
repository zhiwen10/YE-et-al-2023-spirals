function T = getAxonBiasTable(data_folder,save_folder)
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
sensory_label = {'VIS','RSP','AUD','TEa','VISa','VISC','VISrl',...
    'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un','SSs'};
sensory_indx = getRegionIndex(sensory_label,st);                           % get all index values in the cortical sensory area

AUD_label = {'AUD','TEa','VISC'};
AUD_indx = getRegionIndex(AUD_label,st);       % get all index values for auditory area

SSs_label = {'SSs'};
SSs_indx = getRegionIndex(SSs_label,st);                                   % get all index values for SSs area

SSp_label = {'VISrl','SSp-tr','SSp-ll','SSp-ul','SSp-m',...
    'SSp-n','SSp-bfd','SSp-un'};                                           
SSp_indx = getRegionIndex(SSp_label,st);                                   % get all index values for SSp area
% define each single region with axon arbor exclusion criteria
axon_indx{3} = AUD_indx;
axon_indx{4} = AUD_indx;
for kk = [1,2,5,6,7]
    axon_indx{kk} = sensory_indx;
end
for kk = 8:14
    axon_indx{kk} = SSp_indx;
end
axon_indx{15} = SSs_indx;
%% get axon terminals for each cell, grouped by sensory regions
hemi = 1;                                                                  % only look at left hemisphere
for k = 1:numel(sensory_label)
    iregion = sensory_label(k);
    st_region_indx = getRegionIndex(iregion,st);                           % get brain region index
    [cell_id] = getCellPos(allCoords,av,st_region_indx);                   % get soma position in 3d atlas space
    [soma_all2{k},axon_terminal_all2{k}] = getAxonTerminal(...
        cell_id,allCoords,av,axon_indx{k}, hemi);                          % get axon terminals in specifed region boundries
    labels{k} = repmat(string(iregion),numel(soma_all2{k}),1);             % give an area label for each soma
end
%% concatenate all cells
soma_all3 = cat(2, soma_all2{:});
axon_terminal_all3 = cat(2,axon_terminal_all2{:});
labels_all = cat(1,labels{:});
a = cellfun(@(x) not(isempty(x)), axon_terminal_all3);
labels_all = labels_all(a);
soma_all3 = soma_all3(a);
axon_terminal_all3 = axon_terminal_all3(a);
%%
soma_all3a = cellfun(@(x) x(1,:), soma_all3,'UniformOutput',false);
soma_all4 = cat(1,soma_all3a{:});
%%
names = strings(size(soma_all4,1),1);
for i = 1:size(soma_all4,1)
    indxa(i,1) = av(soma_all4(i,1),soma_all4(i,2),soma_all4(i,3));
    names(i,1) = st{st.index == indxa(i,1)-1,4};
end
%%
pat = ["1","2/3","4","5","6"];
for i = 1:5
    pat1 = pat(i);
    TF(:,i) = contains(names,pat1);
end
%% main algorithm for axon terminal svd
[soma_center,axon_vector,angle1,polarity] = getAxonSVD(...
    axon_terminal_all3,soma_all3);   
%% calculate bias angle for each cell, reference to SSp center
clear ang_diff1 ang_diff2
center = [244,542]; % SSp-un
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
    angle1,polarity,ang_diff2,labels_all, TF);
%%
T = renamevars(T,["z","axon_vector","angle1","polarity","ang_diff2","TF"], ...
                 ["soma_polar_angle","axon_bias","bias_angle",...
                 "pc_ratio","center_bias_angle","layers"]);
%%
layer_counts = sum(T{:,9},1);
%%
writetable(T,fullfile(save_folder,'Axon_bias_all_cells3.csv'));