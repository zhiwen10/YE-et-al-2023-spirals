
%% navigating to swc files of interest
global githubDir
folder = 'E:\BIL\download.brainimagelibrary.org+8811\0f\cd\0fcde5fdd6f7ccb2';
cd(folder)

dataDir = dir(folder);
for ii = 3:numel(dataDir)
    swcFileDir{ii-2} = dir(fullfile(folder, dataDir(ii).name, '*reg.swc')); % only look for files registered to Allen CCF
end


%% loading coords of each neuron from each mouse
ctr=0;
for mn = 1:numel(swcFileDir)
    for ii = 1:numel(swcFileDir{mn})
        ctr = ctr+1;
        thisFile = fullfile(swcFileDir{mn}(ii).folder, swcFileDir{mn}(ii).name);
        T = readtable(thisFile,'FileType','text');
        allCoords{ctr} = uint16(round([T.type T.x T.y T.z])); % only keeping x y z coords and type of neuron 
        clear T thisFile
    end
end

%% loading annotation volume and structure tree
githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'allenCCF')))
av = readNPY('annotation_volume_10um_by_index.npy');
st = loadStructureTree();

%% finding neurons with axons that pass through a structure of interest

% regions of interest where axon terminates or passes through
axonRegions = {"MRN", "SCm"}; 
axonRegInd = [];
for r = 1:numel(axonRegions)
    axonRegID = st.id(strcmp(st.acronym, axonRegions{r}));
    axonRegInd = [axonRegInd; find(cellfun(@(x)contains(x, sprintf('/%d/', axonRegID)), st.structure_id_path))];
end

% regions of interest where soma lies
somaRegions =  {"MOp", "MOs"}; % {"SSp"};
somaRegInd = [];
for r = 1:numel(somaRegions)
    somaRegID = st.id(strcmp(st.acronym, somaRegions{r}));
    somaRegInd = [somaRegInd; find(cellfun(@(x)contains(x, sprintf('/%d/', somaRegID)), st.structure_id_path))];
end

ctr=0; 
clear projAxonsInROI projectingSomaLocations projectingNeuronIDs
for nIdx = 1:numel(allCoords)
    somaCoords = allCoords{nIdx}(:,1)==1;
    axonCoords = allCoords{nIdx}(:,1)==2;
    coordsIdx = round(allCoords{nIdx}(:,2:4)/10); % compatibility with av (10 micron resolution)
    % fixing rounding errors:
    coordsIdx(coordsIdx<1) = 1;
    coordsIdx(coordsIdx(:,1)>size(av,1),1) = size(av,1); coordsIdx(coordsIdx(:,2)>size(av,2),2) = size(av,2); coordsIdx(coordsIdx(:,3)>size(av,3),3) = size(av,3); 
    avIdxAxon = sub2ind(size(av), coordsIdx(axonCoords,1),  coordsIdx(axonCoords,2),  coordsIdx(axonCoords,3));
    stIdxAxon = av(avIdxAxon);
    avIdxSoma = sub2ind(size(av), coordsIdx(somaCoords,1),  coordsIdx(somaCoords,2),  coordsIdx(somaCoords,3));
    stIdxSoma = av(avIdxSoma);
    % seeing whether roi ever contains axon and keeping axon coords that pass through ROI
    if any(ismember(stIdxAxon, axonRegInd)) && any(ismember(stIdxSoma, somaRegInd))
        ctr = ctr+1;
        projectingNeuronIDs(ctr) = nIdx;
        projectingSomaLocations(ctr) = unique(st.name(stIdxSoma));
        [projAxonsInROI{ctr}(:,1), projAxonsInROI{ctr}(:,2), projAxonsInROI{ctr}(:,3)] = ind2sub(size(av), avIdxAxon(ismember(stIdxAxon, axonRegInd)));
    end
    clear axonCoords coordsIdx stIdxAxon avIdxAxon stIdxSoma avIdxSoma
end

% allNames = unique(st.name(stIdx))

%% plotting projecting neurons of interest

% coloring somas as a gradient across A/P, M/L, or D/V axes
clear apCoord mlCoord dvCoord
for ii = 1:numel(projectingNeuronIDs)
    nIdx = projectingNeuronIDs(ii);
    type = 1; % interested in somas only
    theseCoords = allCoords{nIdx}(:,1)==type;
    apCoord(ii) = min(allCoords{nIdx}(theseCoords,2));
    mlCoord(ii) = min(allCoords{nIdx}(theseCoords,4));
    dvCoord(ii) = min(allCoords{nIdx}(theseCoords,3));
    clear theseCoords
end

% to color somas and their corresponding axons by position of soma along
% ap, ml, or dv axes; uncomment the 3 lines of code for the axis you want
% to look at
clear colorIdx
apCoord = double(apCoord);
colorIdx = round((apCoord-min(apCoord))/(max(apCoord)-min(apCoord))*256); % standardizing distribution of soma ap coords to 0-256 range
colorIdx(colorIdx==0)=1;

% mlCoord = double(mlCoord);
% colorIdx = round((mlCoord-min(mlCoord))/(max(mlCoord)-min(mlCoord))*256); % standardizing distribution of soma ml coords to 0-256 range
% colorIdx(colorIdx==0)=1;
 
% dvCoord = double(dvCoord);
% colorIdx = round((dvCoord-min(dvCoord))/(max(dvCoord)-min(dvCoord))*256); % standardizing distribution of soma dv coords to 0-256 range
% colorIdx(colorIdx==0)=1;

plotBrainGrid; hold on;
colors = colorcet('R2');
for ii = 1:numel(projectingNeuronIDs)
    nIdx = projectingNeuronIDs(ii);
    somaCoords = allCoords{nIdx}(:,1)==1;
    axonCoords = allCoords{nIdx}(:,1)==2;
    plot3(allCoords{nIdx}(somaCoords,2)/10, allCoords{nIdx}(somaCoords,4)/10, allCoords{nIdx}(somaCoords,3)/10, '.', 'Color', colors(colorIdx(ii),:), 'MarkerSize', 12); % y and z coords need to be flipped around here, I think this may have to do with plot3 and BrainGrid compatibility?
    plot3(projAxonsInROI{ii}(:,1), projAxonsInROI{ii}(:,3), projAxonsInROI{ii}(:,2), '.', 'Color', colors(colorIdx(ii),:), 'MarkerSize', 5); 
    hold on;
end
axis equal
% set(gca,'Xdir','reverse')


%% plotting a single neuron and comparing to 3D morphology in paper as a sanity check
% can compare plot here to corresponding neurons found here: https://braintell.org/projects/fullmorpho/

% f = figure; 
noi = allCoords{349}; % neuron of interest: brain ID 210254, neuron ID *20600
plotBrainGrid; hold on;
colors = get(gca, 'ColorOrder');
for type = 1:max(noi(:,1))
    theseCoords = noi(:,1)==type;
    
    plot3(noi(theseCoords,2)/10, noi(theseCoords,4)/10, noi(theseCoords,3)/10, '.', 'Color', colors(type,:)); % y and z coords need to be flipped around here, I think this may have to do with plot3 and BrainGrid compatibility?
    hold on;
end
axis equal
