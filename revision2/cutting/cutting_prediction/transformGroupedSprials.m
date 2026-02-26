function spirals = transformGroupedSprials(archiveCell,tform)
%% raw spirals 
spiral_length = cellfun(@(x) size(x,1), archiveCell);
spiral_sequence = archiveCell(spiral_length>=2);
spirals= cell2mat(spiral_sequence);
spirals(spirals(:,4)==0,:) = -1;
[spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
    tform,spirals(:,1),spirals(:,2));
spirals(:,1:2) = round(spiralsT); 