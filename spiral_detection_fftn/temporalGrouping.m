%% temporal grouping
pwAll =unique(pwAll, 'rows');
pwAll = sortrows(pwAll,5);
%%
allFrames = unique(pwAll(:,end));
firstFrame = allFrames(1);
lastFrame = allFrames(end);
frameIterator =  lastFrame -firstFrame+1;
%%
% indx1 = find(allFrames == 65009);
indx1  = 1;
cell1 = {};
cFrame = allFrames(indx1);
%%
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
        currentSpirals = pwAll(pwAll(:,end)==cFrame,:);
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
    nextSpirals = pwAll(pwAll(:,end)==nextFrame,:);
    if not(isempty(nextSpirals)) & not(isempty(cell1))
        [archiveCell,cell1,currentSpirals,newRound,newAdd] = groupSpirals(archiveCell,cell1,nextSpirals,newRound,newAdd);
        % currentSpirals = currentSpirals;
    else
        currentSpirals = [];
    end   
    cFrame = cFrame+1;  
    T1(i) = toc;
    if mod(i,1000)==0
        fprintf('Frame %g / totalFrame %g ; time elapsed %g seconds \n', [i,frameIterator,T1(i)-T1(i-999)]);
    end   
end
%% expected spirals numbers
% [pwIndx] = ismember(pwAll(:,end),allFrames(indx1):allFrames(indx1)+numel(allFrames)-1);
[pwIndx] = ismember(pwAll(:,end),allFrames(indx1):allFrames(indx1)+frameIterator);
epochSpirals = pwAll(pwIndx,:);
% grouped + ungrouped spiral numbers
groupedSpirals = [cell2mat([archiveCell;cell1]);currentSpirals];
groupedSpirals = sortrows(groupedSpirals,4);
%% check duplicates
[u,I,J] = unique(groupedSpirals, 'rows', 'first');
ixDupRows = setdiff(1:size(groupedSpirals,1), I);
dupRowValues = groupedSpirals(ixDupRows,:);
%% check missed spirals
[a,b] = ismember(pwAll,groupedSpirals,'rows');
missed = pwAll(~a,:);
%%
indx4 = cellfun(@(x) size(x,1), archiveCell);
figure; 
histogram(indx4)
%%
% figure; 
% for i = 1:150
%     subplot(15,10,i)
%     pwAllEpoch1 = archiveCell{i};
%     scatter(pwAllEpoch1(:,1),pwAllEpoch1(:,2),[],pwAllEpoch1(:,4),'filled');   
%     set(gca,'Ydir','reverse');
%     % set(gca, 'XTickLabel', '', 'YTickLabel', '')
%     xlim([0 512]); ylim([0 512]);
%     % axis equal;
%     % axis off; 
%     % box off;
% end
%%
indx2 = cellfun(@(x) size(x,1)>3, archiveCell);
groupedCells = archiveCell(indx2);
%%
filteredSpirals = cell2mat(groupedCells);
%%
figure; histogram(filteredSpirals(:,3));
% ylim([0,30000]);
%% kernel density plot
figure; 
ax2 = subplot(1,1,1);
% imagesc(mimg);
% hold on
scatter_kde(filteredSpirals(:,1),filteredSpirals(:,2),'filled', 'MarkerSize', 5')
set(gca,'Ydir','reverse')
set(ax2, 'XTickLabel', '', 'YTickLabel', '');

xlim([0 512]); ylim([0 512]);
axis off; axis image
box off
colormap(ax2,parula)
%%
meanSpirals = cellfun(@(x) mean(x,1), groupedCells,'UniformOutput',false);
meanSpirals = cell2mat(meanSpirals);
indx5 = cellfun(@(x) size(x,1), groupedCells);
%% kernel density plot
figure; 
ax2 = subplot(1,1,1);
% imagesc(mimg);
% hold on
scatter_kde(meanSpirals (:,1),meanSpirals (:,2),'filled', 'MarkerSize', 5')
set(gca,'Ydir','reverse')
set(ax2, 'XTickLabel', '', 'YTickLabel', '');
xlim([0 512]); ylim([0 512]);
axis off; axis image
box off
colormap(ax2,parula)
%% scatter size plot
indx5 = log10(indx5-min(indx5)+1);
figure; 
ax2 = subplot(1,1,1);
scatter(meanSpirals (:,1),meanSpirals (:,2),3,meanSpirals (:,3),'filled')
% scatter(meanSpirals(:,1),meanSpirals(:,2),3,indx5,'filled')
set(gca,'Ydir','reverse')
set(ax2, 'XTickLabel', '', 'YTickLabel', '');
xlim([0 512]); ylim([0 512]);
axis off; axis image
box off
colormap(ax2,parula)
%%
roi1 = drawpolygon;
%%
tf = inROI(roi1,meanSpirals(:,1),meanSpirals(:,2));
selectedSpirals = groupedCells(tf);
%%
% selectedSpirals = groupedCells(:);
%%
indx3 = cellfun(@(x) size(x,1)>15, selectedSpirals);
selectedSpirals3 = selectedSpirals(indx3);
%%
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat')
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
%%
T = readtable('C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spiralSessions.xlsx');
kk = 21;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
serverRoot = expPath(mn, td, en);
reg_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
dfolder = fullfile(reg_folder,mn,td,num2str(en));
tformName = dir(fullfile(dfolder,'*tform.mat')).name;
load(fullfile(dfolder,tformName));
%%    
scale = 1;
figure;
for i = 1:numel(selectedSpirals3)
    subplot(6,6,i)
    clear pwAllEpoch1 u v
    pwAllEpoch1 = selectedSpirals3{i};
    [u,v] = transformPointsForward(tform,pwAllEpoch1(:,1),pwAllEpoch1(:,2));
    overlayOutlines(coords,scale,[0.8,0.8,0.8]);
    set(gca, 'YDir','reverse');
    scatter(u,v,5,pwAllEpoch1(:,5)-pwAllEpoch1(1,5),'filled');   
    set(gca,'Ydir','reverse');
    % set(gca, 'XTickLabel', '', 'YTickLabel', '')
    % xlim([0 512]); ylim([0 512]);
    axis image; axis off;
    title(['epoch ' num2str(i)]); 
end
%% 
% epoch = [1,2,3,4,6,16,19,22,23,25,26,29,30,31,32,33,34,35];
% epoch = [1,9,11,14,20,26];
epoch = 1:2;
good_epoch = selectedSpirals3(epoch);
epoch_start_end = cellfun(@(x) [x(1,5),x(end,5)], good_epoch,'UniformOutput',false);
epoch_start_end2 = cell2mat(epoch_start_end);
%%
h = figure;
for i = 1:numel(good_epoch)
    subplot(4,5,i)
    clear pwAllEpoch1 u v
    pwAllEpoch1 = good_epoch{i};
    [u,v] = transformPointsForward(tform,pwAllEpoch1(:,1),pwAllEpoch1(:,2));
    overlayOutlines(coords,scale,[0.8,0.8,0.8]);
    set(gca, 'YDir','reverse');
    scatter(u,v,5,pwAllEpoch1(:,5)-pwAllEpoch1(1,5),'filled');   
    set(gca,'Ydir','reverse');
    % set(gca, 'XTickLabel', '', 'YTickLabel', '')
    % xlim([0 512]); ylim([0 512]);
    axis image; axis off;
    title([num2str(epoch_start_end2(i,1)) '-' num2str(epoch_start_end2(i,2))],'fontSize',8); 
end
%%
print(h, 'grouped_epoch', '-dpdf', '-bestfit', '-painters');
%%
% good: 5,8,11,13,15
%%
[U,V,t,mimg] = get_wf_svd1(serverRoot);
dV = [zeros(size(V,1),1) diff(V,[],2)];

%%
for i = 1:numel(epoch)
% for i = 1
    current_epoch = selectedSpirals3{epoch(i)};
    frameStart1 = current_epoch(1,5); frameEnd1 = current_epoch(end,5);
    frameStart = frameStart1-35; frameEnd = frameEnd1+35;
    video_name = [mn '_' tdb '_' num2str(en) '_' num2str(frameStart1) '-' num2str(frameEnd1)];
    make_raw_phase_video_by_frame2(U(:,:,1:50),dV(1:50,:),t,BW1,mimg,projectedAtlas1,coords,tform,frameStart,frameEnd,pwAll,video_name);
    close all
end
%%
U1 = imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
BW = logical(projectedAtlas1(1:4:end,1:4:end));
params.lowpass = 0;
params.gsmooth = 0;
for i = 1:numel(epoch)
% for i = 1
    current_epoch = selectedSpirals3{epoch(i)};
    frameStart = current_epoch(1,5); frameEnd = current_epoch(end,5);
    frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
    dV1 = dV(1:50,frameTemp);
    t1 = t(frameTemp);
    rate1 = 0.1;
    [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U1(1:4:end,1:4:end,1:50),dV1,t1,params,rate1);
    tracePhase1 = tracePhase1(:,:,1+35/rate1:end-35/rate1); % reduce 2*35 frames after filter data 
    tracePhase1 = permute(tracePhase1,[3,1,2]);
    %% phase map video
    fname = [mn '_' tdb '_' num2str(en) '_' num2str(frameStart) '-' num2str(frameEnd)];
    cmap = 'hsv';
    caxis1 = [min(tracePhase1(:)),max(tracePhase1(:))];
    make_widefield_video(tracePhase1,BW,cmap,caxis1,fname);
end
%%
frameID2 = 197:213;
rate1 = 0.1;
[tracePhase, pwAllEpoch, cells] = upsampleEpoch(U1,dV,t,frameID2,params,rate1);
% pwAllEpoch1 = pwAllEpoch(pwAllEpoch(:,3)>50,:);
frameID = 1:size(tracePhase,3);
% save('frame197to213upsampled.mat','frameID','-append')
%%
frameID2 = 6804:6819;
rate1 = 0.1;
[tracePhase, pwAllEpoch, cells] = upsampleEpoch(U1,dV,t,frameID2,params,rate1);
% pwAllEpoch1 = pwAllEpoch(pwAllEpoch(:,3)>50,:);
frameID = 1:size(tracePhase,3);
% save('frame197to213upsampled.mat','frameID','-append')
%%
