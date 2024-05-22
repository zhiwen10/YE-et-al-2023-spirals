function [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals)
%% spatiotemporal grouping sprials into clusters 
% spirals in nearby frames with distance <30 pixels are grouped in a cluster
% if no spirals are found nearby the last spiral in a given cluster for at 
% least 3 frames, then that cluster is saved as an archived cluster cell.

% cluster   a    b    c    d   e     
% frame 1:  1    1    1    0   0    
% frame 2:  0    0    1    1   0     
% frame 3:  0    0    1    1   0       
% frame 4:  0    0    0    1   0     
% frame 5:  0    0    0    0   0     
% frame 6:  0    0    0    0   1      
%% iterate frame by frame
% frame 1: add 3 new spirals in 3 new cell clusters {a1,b1,c1} within the 
% imtermediate cluster assembly <cell1>

% frame 2: if a spiral (c2) in the new frame2 was nearby the last spiral (c1) 
% in the imtermediate cluster <cell1>, then attach that sprial(c2) to the 
% cluster,{a1, b1,(c1,c2)};

% if a spiral (d1) in the new frame2 did not belong to any exisiting cluster,
% create a new cluster, {a1,b1,(c1,c2),d1};
% keep a tally of how many frames each new cluster went through 
% <newRound [2,2,2,1]>, and how many new spirals were added to each cluster 
% <newAdd [1,1,2,1]>.  

% repeat the same process for frame3,4, which created intermediate
% clsuter assembly <cell1> with {a1,b1,(c1,c2,c3),{d2,d3,d4)}, and each
% cluster went through <newRound [4,4,4,3]> frame iterations, and <newAdd
% [1,1,3,3]> sprials were added for each cell cluster. 

% at frame 4, for cluster a, b, at least 3 frames went by without adding new 
% spirals to their cluster, thereby cluster a,b were archived, and dropped 
% from cluster assemby <cell1 {(c1,c2,c3),(d2,d3,d4)}>, <newRound [4,3]>, 
% <newAdd [3,3]>

% interate through all frames to cluster all spirals.
%%
allFrames = unique(filteredSpirals(:,end));                            % get unique frames with spirals
firstFrame = allFrames(1);
lastFrame = allFrames(end);
frameIterator =  lastFrame -firstFrame+1;
%%
indx1  = 1;
cell1 = {}; archiveCell = {};
cFrame = allFrames(indx1); 
tic
for i = 1:frameIterator  
    %%
    % iterate frame by frame
    % at initial step/frame i = 1, when there is 
    % no cluster in the intermediate cluster assembly <cell1>,
    % store all spirals in the current frame 
    % to new clusters in <cell1> 
    %%
    if isempty(cell1) 
        %%
        % except when <cell1> is empty in the middle 
        % of the iteration process(frame i>1),
        % all previous spirals were archived successfully,
        % then currentSpirals only contain spirals that 
        % were not archived yet. 

        % to save time, when archiveCell has more than 100 clusters,
        % only look for currentSpirals that were not archived 
        % in the recent 100 achived clusters
        %%
        if numel(archiveCell)>100  
            tempSpirals = cell2mat(archiveCell(end-100:end));
        else
            tempSpirals = cell2mat(archiveCell);
        end
        currentSpirals = ...
            filteredSpirals(filteredSpirals(:,end)==cFrame,:);

        if not(isempty(currentSpirals)) 
            %%
            % check if the currentSpiral has been grouped 
            % previously already
            % if not, then iterate through the remaining spirals, 
            % otehrwise skip.
            %%
            if not(isempty(tempSpirals)) 
                a = ismember(currentSpirals,tempSpirals,'rows');
                % only look at sprials that were not archived 
                % already in currentSpirals
                if sum(a)>0
                    currentSpirals = currentSpirals(~a,:);
                end
            end
        end
        if not(isempty(currentSpirals))                                     % initialize cell1, newRound and newAdd
            rowDist = ones(1,size(currentSpirals,1));
            cell1 = mat2cell(currentSpirals,rowDist);
            newRound = zeros(numel(cell1),1);
            newAdd = zeros(numel(cell1),1);
        end
    elseif not(isempty(cell1)) & not(isempty(currentSpirals))  
        %%
        % when there are intermediate cluster assembly 
        % <cell1> existing and some spirals in the current frame 
        % were not attached to existing clusters, 
        % then attach these sprials as new clusters in <cell1>, 
        % and add new tally in <newRound> and <newAdd>
        %%
        rowDist = ones(1,size(currentSpirals,1));
        cell2 = mat2cell(currentSpirals,rowDist);
        cell1 = [cell1;cell2];
        newRound = [newRound;zeros(size(currentSpirals,1),1)];
        newAdd = [newAdd;zeros(size(currentSpirals,1),1)];
        currentSpirals = [];
    end 
    newRound = newRound+1;
    nextFrame = cFrame+1;
    nextSpirals = filteredSpirals(filteredSpirals(:,end)==nextFrame,:);
    %%
    % if spirals in the next frame are nearby the last spiral in an
    % exisitng cluster in <cell1>, then attached that spiral in the
    % cluster; otherwise keep that spiral in the currentSpiral list

    % if no new spirals were added to the existing cluster for 
    % more than 3 frames, 
    % then archive that cluster to <archiveCell>
    %%
    if not(isempty(nextSpirals)) & not(isempty(cell1))
        [archiveCell,cell1,currentSpirals,newRound,newAdd] = ...
            groupSpirals(archiveCell,cell1,nextSpirals,newRound,newAdd);
    else
        currentSpirals = [];
    end   
    cFrame = cFrame+1;  
    T1(i) = toc;
    if mod(i,1000)==0
        fprintf(...
            'Frame %g / totalFrame %g ; time elapsed %g seconds \n',...
            [i,frameIterator,T1(i)-T1(i-999)]);
    end   
end
%% check all spirals were archived, and there were no duplication
test_stats = ones(1,2);
[pwIndx] = ismember(filteredSpirals(:,end),...
    allFrames(indx1):allFrames(indx1)+frameIterator);
epochSpirals = filteredSpirals(pwIndx,:);
groupedSpirals = [cell2mat([archiveCell;cell1]);currentSpirals];
groupedSpirals = sortrows(groupedSpirals,4);
[u,I,J] = unique(groupedSpirals, 'rows', 'first');                         % check duplicates
ixDupRows = setdiff(1:size(groupedSpirals,1), I);
dupRowValues = groupedSpirals(ixDupRows,:);
[a,b] = ismember(filteredSpirals,groupedSpirals,'rows');                   % check missed spirals
missed = filteredSpirals(~a,:);
if not(isempty(dupRowValues))
    test_stats(1,1) = 0;
end
if not(isempty(missed))
    test_stats(1,2) = 0;
end  