function [cell_id] = get_cell_id(allCoords,av,st_region_indx)
cell_id = [];
for nIdx = 1:numel(allCoords)
    somaCoords = allCoords{nIdx}(:,1)==1;
    coordsIdx = round(allCoords{nIdx}(:,2:4)/10); % compatibility with av (10 micron resolution)
    avIdxSoma = sub2ind(size(av), coordsIdx(somaCoords,1),  coordsIdx(somaCoords,2),  coordsIdx(somaCoords,3));
    stIdxSoma = av(avIdxSoma);
    if any(ismember(stIdxSoma, st_region_indx))
        cell_id = [cell_id;nIdx];
    end
end