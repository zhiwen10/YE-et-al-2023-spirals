function  [archiveCell,cell1,remain,newRound,newAdd] = ...
    groupSpirals(archiveCell,cell1,nextSpirals,newRound,newAdd)
%%
% check iteratively for each spiral in the next frame (nextSpirals),
% if it is close to the last spiral of any cluster in the intermedate 
% cell assembly <cell1>

% if it is within 30 pixels distance, then attached it to the existing
% cluster in  <cell1>, otherwise keep it in the remaining spiral list as 
% output and create new sprial cluster with this remaining spiral list 
% outside of this function later

% if no new spirals were added to the existing cluster in <cell1> for at
% least 3 frames, then that cluster is archived in the <archiveCell>
%%
% get the last spiral from each cluster in the imtermediate cell assembly
% <cell1>
grouped = [];
lastSpiralsG = cell2mat(cellfun(@(x) x(end,:),cell1,'UniformOutput',false)); 
% if the last spiral was within 3 frames of the next spiral to be checked,
% then consider it to be grouped if distance condition is met, otherwise
% don't even consider
indx1 = (lastSpiralsG(:,end) < nextSpirals(1,end)-2);
lastSpiralsG(indx1,:)=nan;
if not(all(isnan(lastSpiralsG)))
    for j = 1:size(nextSpirals,1)
        tSpirals = nextSpirals(j,:); 
        if not(isnan(tSpirals))
            distance = vecnorm(lastSpiralsG(:,1:2)-tSpirals(1:2),2,2);     % check Euclidean distance
            [minV,indx] = min(distance,[],'omitnan');
            if minV<30
                cell1{indx} = [cell1{indx};tSpirals];                      % attach that spiral in the next frame to existing cluster
                newAdd(indx) = newAdd(indx)+1;                             % tally new addtion
                grouped = [grouped;j];                                     % keep track of the spirals that were grouped
            end 
        end
    end       
end
remain = nextSpirals(setdiff(1:end,grouped),:);                            % spirals that were not grouped, remains
% after more than 3 new frames iterated, 
% and more than 3 frames without new spirals, 
% then drop out that cluster in <cell1> to archiveCell 
% and drop the tally in <newRound> and <newAdd>
newMiss = newRound-newAdd;
indxG =  (newRound>=3 & newMiss>=3);
if any(indxG)
    archiveCell = [archiveCell;cell1(indxG)];
    cell1 = cell1(~indxG);
    newRound = newRound(~indxG);
    newAdd = newAdd(~indxG);
end
end