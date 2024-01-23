function  [archiveCell,cell1,remain,newRound,newAdd] = groupSpirals(archiveCell,cell1,nextSpirals,newRound,newAdd)
grouped = [];
% only use the spiral in the lastest frame for each group
lastSpiralsG = cell2mat(cellfun(@(x) x(end,:),cell1,'UniformOutput',false));
indx1 = (lastSpiralsG(:,end) < nextSpirals(1,end)-2);
lastSpiralsG(indx1,:)=nan;
if not(all(isnan(lastSpiralsG)))
    for j = 1:size(nextSpirals,1)
        tSpirals = nextSpirals(j,:); 
        if not(isnan(tSpirals))
            distance = vecnorm(lastSpiralsG(:,1:2)-tSpirals(1:2),2,2);
            [minV,indx] = min(distance,[],'omitnan');
            if minV<30
                cell1{indx} = [cell1{indx};tSpirals];
                newAdd(indx) = newAdd(indx)+1;
                grouped = [grouped;j];
            end 
        end
    end       
end
remain = nextSpirals(setdiff(1:end,grouped),:);
% lastSpiral = remain;
%%
newMiss = newRound-newAdd;
indxG =  (newRound>=3 & newMiss>=3);
if any(indxG)
    archiveCell = [archiveCell;cell1(indxG)];
    cell1 = cell1(~indxG);
    newRound = newRound(~indxG);
    newAdd = newAdd(~indxG);
end
end