function getWhiskerSingleTrials(data_folder,save_folder)
%% load atlas brain horizontal projection and outline  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
BW2 = BW(1:2:end,1:2:end);
%%
spiral_folder = fullfile(data_folder,'whisker','spirals_all','spirals_grouping');
T = readtable(fullfile(data_folder,'tables','whisker_stim_all2.xlsx'));
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
%% right SSp index
clear areaPath
spath = string(st.structure_id_path);
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
%%
sensoryArea = [];
Utransformed = projectedAtlas1;
scale = 1;
hemi = 'right';
[indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row,col] = ind2sub(size(projectedAtlas1),indexSSp);
indexSSp_right = [col,row];

hemi = 'left';
[indexSSp2,UselectedSSp2] = select_area(sensoryArea,spath,st,coords,...
    Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
[row2,col2] = ind2sub(size(projectedAtlas1),indexSSp2);
indexSSp_left = [col2,row2];
%%
kk = 4;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);
block_en = T.block_folder(kk);
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
% load data and block
session_root = fullfile(data_folder,'whisker','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
U = U./mimg;
bfd(1,1) = T{kk,7};
bfd(1,2) = T{kk,8};
expDir = dir(fullfile(session_root,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
%
sigName = 'rewardValve';
load(fullfile(session_root,[sigName '_raw.mat']));                     % load pd
load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));     % load tlTimes
tt = tsToT(tlTimes, numel(pd)); 
%
[allPD,flipsUp,flipsDown] = schmittTimes(tt, pd, [0.5 0.8]);
flipsUp_iti = diff(flipsUp);
flipsUp_iti = [0;flipsUp_iti];
flipsUp = flipsUp(flipsUp<2700);
%%
load(fullfile(data_folder,'whisker','rfmap',[fname '.mat']));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = mimgtransformed(1:params.downscale:end,1:params.downscale:end);
%%
for ii = 1:numel(flipsUp)
    stimOn(ii,1) = find(t-flipsUp(ii)>0, 1, 'first');
end
frames1 = [-70:70];
frames  = repmat(frames1,[numel(stimOn),1]);
frames_stimOn = frames+stimOn;   
%%
clear Va indx
for i = 1:numel(flipsUp)
    ta = t-flipsUp(i);
    indx(i,1) = find(ta>0, 1,'first');
    tt1(i,:) = t(indx(i,1)-70:indx(i,1)+70)-flipsUp(i);
    Va(:,:,i) = dV(:,indx(i,1)-70:indx(i,1)+70);
end
Ut= imwarp(U,tform,'OutputView',imref2d(size(projectedTemplate1)));
Ut1 = Ut(1:8:end,1:8:end,:);
va1 = reshape(Va,size(Va,1),size(Va,2)*size(Va,3));
wf = reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3))*va1;
wf = reshape(wf, size(Ut1,1),size(Ut1,2),size(Va,2),numel(flipsUp));  
if strcmp(T.hemisphere{kk},'left')        
    wf = flip(wf,2);
end
%% load spirals
load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));    
spiral_duration = cellfun(@(x) size(x,1), archiveCell);
indx2 = (spiral_duration>=2);
groupedCells = archiveCell(indx2);
pwAll = cell2mat(groupedCells); 
% organize sprials by frame
[pwAll(:,1),pwAll(:,2)] = transformPointsForward(...
    tform,pwAll(:,1),pwAll(:,2)); 
pwAll(:,1:2) = round(pwAll(:,1:2));
[lia,locb] = ismember(pwAll(:,1:2),brain_index,'rows');
pwAll = pwAll(lia,:);
if strcmp(T.hemisphere{kk},'left')        
   pwAll(:,1) = 1140-pwAll(:,1);
   pwAll(:,4) = -pwAll(:,4);
end
[lia2,locb2] = ismember(pwAll(:,1:2),indexSSp_right,'rows');
pwAll2 = pwAll(lia2,:);
spiral_cell2 = cell(size(frames_stimOn));
for j = 1:size(frames_stimOn,1)
    for k = 1:size(frames_stimOn,2)
        clear indx1 current_spirals
        current_frame = frames_stimOn(j,k);
        indx1 = find(pwAll2(:,5) == current_frame);
        current_spirals = pwAll2(indx1,:);
        spiral_cell2{j,k} = current_spirals;
    end
end  
frames_post = 71:85;
spiral_large_count = zeros(size(spiral_cell2,1),1);
for i = 1:size(spiral_cell2,1)
    spiral_temp = cat(1,spiral_cell2{i,frames_post});
    spiral_large_count(i,1) = sum(spiral_temp(:,3)>=60);
end
index = (spiral_large_count>=5);
trial_select = find(index);
%%
% indx2 = [6,21,36,37];
% indx2 = [6,36];
% indx2 = [21,6,31];
% wf1 = squeeze(wf(:,:,:,trial_select(indx2)));
wf1 = squeeze(wf(:,:,:,[81,486,31]));
% wf2 = imresize(wf1,[660,570]);
save(fullfile(save_folder,'single_trial_maps2.mat'),'wf1');