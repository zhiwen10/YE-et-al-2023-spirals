function h2 = plotWfTraceSingleTrial(data_folder,T1,session,trial_indx)
%% load atlas brain horizontal projection and outline 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
downscale = 16;
BW3 = BW(1:downscale:end,1:downscale:end);
%%
pixel(1,:) = [870,850]; % VISp
pixel(2,:) = [775,650]; % RSP
pixel(3,:) = [590,750]; % SSp-ul
pixel(4,:) = [520,850]; % SSp-ll
pixel(5,:) = [480,960]; % SSp-m
pixel(6,:) = [550,960]; % SSp-n
pixel(7,:) = [660,950]; % SSp-bfd
pixel = round(pixel/downscale);
%% session info
mn = char(T1.MouseID{session});
tda = T1.date(session);  
en = T1.folder(session); 
block_en = T1.block_folder(session); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%%
% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
expDir = dir(fullfile(session_root,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
% get photodiode time
win = [0,3000];
allPD2 = getPhotodiodeTime(session_root,win);
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1)';
allPD2 = allPD2(1:ntrial);
allPD2 = allPD2(norepeat_indx);
%% register to atlas
load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(dV(1:50,:));
%% get time index for all events
for i = 1:numel(allPD2)
    ta  = t-allPD2(i);
    indx2(i,1) = find(ta>0, 1,'first');
    t2(i,1) = t(indx2(i,1));
end
%% plot
indx_t = indx2(trial_indx);                                                      % stimOn time index in time stamp
Ut1a = double(reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3)));
trace = Ut1a*dV1(:,indx_t-17:indx_t+35);
trace2d = reshape(trace,size(Ut1,1),size(Ut1,2),size(trace,2))./mimgt;
%%
color2 = cbrewer2('seq','YlOrRd',9);
nameList = {'VISp','RSP','SSp_ul','SSp_ll','SSp_m','SSp_n','SSp_bfd'};
tstep = [-0.5:1/35:1];
h2 = figure('Renderer', 'painters', 'Position', [50 50 300 300]);
ax2 = subplot(1,1,1);
for i = 1:7
    trace_temp =squeeze(trace2d(pixel(i,1),pixel(i,2),:));
    plot(tstep,trace_temp+0.02*i,'color',color2(i+2,:));
    hold on;
end    
set(ax2,'YTickLabel',[]);
for i = 1:7
    text(-0.7,0.02*i,nameList{i},'Color', color2(i+2,:),'Interpreter','None');
end
xlabel('Time (s)','fontSize',9);
line_ax1 = xline(0,'--');