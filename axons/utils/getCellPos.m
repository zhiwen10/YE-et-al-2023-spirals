function [cell_id] = getCellPos(allCoords,av,filter_region)
cell_id = [];
for nIdx = 1:numel(allCoords)
    somaCoords = allCoords{nIdx}(:,1)==1;                                  % soma positions are [type1] in [column1]
    coordsIdx = round(allCoords{nIdx}(:,2:4)/10);                          % downsample 3d coordinates to 10um resolution
    avIdxSoma = sub2ind(size(av), coordsIdx(somaCoords,1),  ...
        coordsIdx(somaCoords,2),  coordsIdx(somaCoords,3));                % find out the soma index in 3d atlas space
    stIdxSoma = av(avIdxSoma);                                             % find out the soma key-value in 3d atlas space 
    if any(ismember(stIdxSoma, filter_region))
        cell_id = [cell_id;nIdx];                                          % if soma coordiates are in the specified region, then save
    end
end