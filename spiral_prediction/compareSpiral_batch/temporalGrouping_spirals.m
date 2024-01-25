% add path for dependencies 
githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'spikes')))
addpath(genpath(fullfile(githubdir2, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubdir2, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubdir2, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
addpath(genpath(fullfile(githubdir2, 'npy-matlab')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow1';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T1 = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
%%
dfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_raw';
savefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_raw_fftn';
% kk = 16; % zye41 STR
% kk = 11; % zye60
% kk = 39; % zye67
test_stats(size(T1,1),1) = 0;
% for kk = 1:size(T1,1) % zye58
for kk = 28
    ops = get_session_info(T1,kk);
    % kk = 11;
    clear spiralsT
    mn = T1.mouseId{kk};
    tda = T1.date(kk);
    en = T1.en(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    serverRoot = expPath(mn, td, en);
    %%
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(dfolder,[fname '_spirals.mat']));
    filteredSpirals = pwAll(pwAll(:,3)>=40,:);
    %% temporal grouping
    filteredSpirals =unique(filteredSpirals, 'rows');
    filteredSpirals = sortrows(filteredSpirals,5);
    %%
    allFrames = unique(filteredSpirals(:,end));
    firstFrame = allFrames(1);
    lastFrame = allFrames(end);
    frameIterator =  lastFrame -firstFrame+1;
    %%
    indx1  = 1;
    cell1 = {};
    cFrame = allFrames(indx1);
    archiveCell = {};
    tic
    for i = 1:frameIterator 
    % for i = 1:100      
    %     if cFrame == 9640
    %         fprintf('pause');
    %     end
        if isempty(cell1) 
            if numel(archiveCell)>100
                tempSpirals = cell2mat(archiveCell(end-100:end));
            else
                tempSpirals = cell2mat(archiveCell);
            end
            currentSpirals = filteredSpirals(filteredSpirals(:,end)==cFrame,:);
            if not(isempty(currentSpirals)) 
                % check if the currentSpiral has been grouped previously already
                % if not, then iterate through the remaining spirals, otehrwise skip.
                if not(isempty(tempSpirals)) 
                a = ismember(currentSpirals,tempSpirals,'rows');
                    if sum(a)>0
                        currentSpirals = currentSpirals(~a,:);
                    end
                end
    %             if not(isempty(nextSpirals)) & isempty(remain) 
    %                 currentSpirals = [];
    %             end
            end
            if not(isempty(currentSpirals)) 
                % initialize cell1, newRound and newAdd
                rowDist = ones(1,size(currentSpirals,1));
                cell1 = mat2cell(currentSpirals,rowDist);
                newRound = zeros(numel(cell1),1);
                newAdd = zeros(numel(cell1),1);
            end
        elseif not(isempty(cell1)) & not(isempty(currentSpirals))  
            rowDist = ones(1,size(currentSpirals,1));
            cell2 = mat2cell(currentSpirals,rowDist);
            cell1 = [cell1;cell2];
            newRound = [newRound;zeros(size(currentSpirals,1),1)];
            newAdd = [newAdd;zeros(size(currentSpirals,1),1)];
            currentSpirals = [];
        end 
        % fprintf('Cells %g; newRound %g; newAdd %g\n', [numel(cell1),numel(newRound),numel(newAdd)]);
        newRound = newRound+1;
        %%
        nextFrame = cFrame+1;
        nextSpirals = filteredSpirals(filteredSpirals(:,end)==nextFrame,:);
        if not(isempty(nextSpirals)) & not(isempty(cell1))
            [archiveCell,cell1,currentSpirals,newRound,newAdd] = groupSpirals(archiveCell,cell1,nextSpirals,newRound,newAdd);
            % currentSpirals = currentSpirals;
        else
            currentSpirals = [];
        end   
        cFrame = cFrame+1;  
        t1(i) = toc;
        if mod(i,1000)==0
            fprintf('Frame %g / totalFrame %g ; time elapsed %g seconds \n', [i,frameIterator,t1(i)-t1(i-999)]);
        end   
    end
    %%
    savefile = fullfile(savefolder,[fname '_spirals_group_fftn.mat']);    
    save(savefile,'archiveCell');
    %% expected spirals numbers
    [pwIndx] = ismember(filteredSpirals(:,end),allFrames(indx1):allFrames(indx1)+frameIterator);
    epochSpirals = filteredSpirals(pwIndx,:);
    groupedSpirals = [cell2mat([archiveCell;cell1]);currentSpirals];
    groupedSpirals = sortrows(groupedSpirals,4);
    % check duplicates
    [u,I,J] = unique(groupedSpirals, 'rows', 'first');
    ixDupRows = setdiff(1:size(groupedSpirals,1), I);
    dupRowValues = groupedSpirals(ixDupRows,:);
    % check missed spirals
    [a,b] = ismember(filteredSpirals,groupedSpirals,'rows');
    missed = filteredSpirals(~a,:);
    if not(isempty(dupRowValues))
        test_stats(kk,1) = 1;
    end
    if not(isempty(missed))
        test_stats(kk,1) = 1;
    end    
end

