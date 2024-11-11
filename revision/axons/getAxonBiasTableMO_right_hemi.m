function getAxonBiasTableMO_right_hemi(data_folder,save_folder)
%% get axon arbor bias for each neuron in sensory region with SVD
% first, define local axon regions, so that long range projections are excluded;
% then, do svd on axon terminal positions for each neuron.

% if soma is in VIS or RSP, no restriction;
% if soma is in auditory cortex, or secondary SS, restrict to within the
% area, because of cross-modality projections;
% if soma is in somatosensory cortex, restrict to the large SSp region, 
% so that long-range motor projections are excluded.
%% load singel cell morphology dataset
load(fullfile(data_folder,'axons\all_cell_with_parents.mat'));
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
av = readNPY(fullfile(data_folder,'tables',...
    'annotation_volume_10um_by_index.npy'));
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%%
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
%% define axon arbor search regions by its soma location
MO_label = {'MOp','MOs','ACA','PL','ILA','ORB'};
MO_indx = getRegionIndex(MO_label,st);                           % get all index values in the cortical sensory area

MOp_label = {'MOp'};
MOp_indx = getRegionIndex(MOp_label,st);       % get all index values for auditory area

MOs_label = {'MOs'};
MOs_indx = getRegionIndex(MOs_label,st);                                   % get all index values for SSs area

ACA_label = {'ACA','PL','ILA','ORB'};                                           
ACA_indx = getRegionIndex(ACA_label,st);                                   % get all index values for SSp area
% define each single region with axon arbor exclusion criteria
% axon_indx{1} = MOp_indx;
% axon_indx{2} = MOs_indx;
% for k = 3:6
%     axon_indx{k} = ACA_indx;
% end
for k = 16:21
    axon_indx{k} = MO_indx;
end
%%
all_label = {'VIS','RSP','AUD','TEa','VISa','VISC','VISrl',...
    'SSp-tr','SSp-ll','SSp-ul','SSp-m','SSp-n','SSp-bfd','SSp-un','SSs',...
    'MOp','MOs','ACA','PL','ILA','ORB'};
%% get axon terminals for each cell, grouped by sensory regions
hemi = 2;                                                                  % only look at left hemisphere
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
%%
% theta = 90;
% for i = 1:size(axon_vector,1)
%     v1 = axon_vector(i,:);
%     vr = v1*[1;1i]*exp(-1i*theta*pi/180);
%     angle1(i,1) = round(rad2deg(atan(imag(vr)/real(vr))))+1;
% end
% angle1(angle1>90) = 90; % one cell came out at 91, let's round it to 90 for color interpolation;
% figure;
% histogram(angle1,20);
%% calculate bias angle for each cell, reference to SSp center
% center = [244,542]; % SSp-un
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
writetable(T,fullfile(save_folder,'Axon_bias_all_cells_MO.csv'));