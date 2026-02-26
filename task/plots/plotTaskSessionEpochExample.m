function hs14b = plotTaskSessionEpochExample(data_folder,save_folder)
%%
win = [0,5000];
downscale = 16;
Visp(1,:) = [810,280]; % left hemisphere
Visp(2,:) = [810,880]; % right hemisphere
Visp(3,:) = [900,280]; % left hemisphere
Visp(4,:) = [900,880]; % right hemisphere
Visp(5,:) = [512,850]; % right SSp-ul
Visp(6,:) = [768,640]; % right RSP
Visp_ds = round(Visp./downscale);
%%
mn = 'ZYE_0091';
T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
T1 = T_session(T_session.label == "task",:);
T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
kk = 1;
load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome.mat']));
load(fullfile(data_folder,'task','task_outcome',[mn '_task_trial_ID.mat']));
indexa = (trial_all.session == kk);
T_all = T_all(indexa,:);
mn = char(T1.MouseID{kk});
tda = T1.date(kk);  
en = T1.folder(kk); 
block_en = T1.block_folder(kk); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task','task_svd',fname);
[U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
expDir = dir(fullfile(session_root,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
% get photodiode time
% serverRoot = expPath(mn,td,en);
% win = [0,3000];
allPD2 = getPhotodiodeTime(session_root,win);
%% rotaryEncoder
sigName = 'rotaryEncoder';
load(fullfile(session_root,[sigName '_raw.mat']));                      % load pd
load(fullfile(session_root,[sigName '_timestamps_Timeline.mat']));         % load tlTimes
% tlFile = fullfile(session_root, [sigName '.raw.npy']); 
% pd = readNPY(tlFile);
% tlFile = fullfile(session_root, [sigName '.timestamps_Timeline.npy']);
% tlTimes = readNPY(tlFile);
tt = tsToT(tlTimes, numel(pd));  
fs1 = 1/mean(diff(tt));
wh = correctCounterDiscont(pd);
vel = computeVelocity2(wh, 0.01, fs1); 
vela = interp1(tt,vel,t); 
% vel_abs = abs(vela);
% vel_abs = vel_abs-mean(vel_abs);
% vel_abs = vel_abs';
%%
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1)';
response_time = block.events.responseTimes(1:ntrial)-block.events.stimulusOnTimes(1:ntrial);
%%
allPD2 = allPD2(1:ntrial);
left_contrast = block.events.contrastLeftValues(1:ntrial)';
right_contrast = block.events.contrastRightValues(1:ntrial)';
contrast_all = [left_contrast, right_contrast];
contrasts_unique = contrast_all(norepeat_indx,:);
% load atlas registration
load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(V(1:50,:));
% dV1 = double(dV(1:50,:));
%% get all traces for 4 pixels
wf1 = [];
for j = 1:6
    wfa = squeeze(Ut1(Visp_ds(j,1),Visp_ds(j,2),:))'*dV1;
    wfa = wfa./mimgt(Visp_ds(j,1),Visp_ds(j,2));
    wf1 = [wf1;wfa];
end
wf1 = wf1*100;
%
wf2 = wf1';
trace1_demean = wf2-mean(wf2,1);
trace1_demean = double(trace1_demean);
%% get time index for all events
clear indx
for i = 1:numel(allPD2)
    ta  = t-allPD2(i);
    indx(i,1) = find(ta>0, 1,'first');
    t2(i,1) = t(indx(i,1));
end
%%
contrast_all1 = contrast_all;
contrast_all1(:,1) = contrast_all1(:,1)*-1;
contrast_all1 = sum(contrast_all1,2);
contrast_ttl1 = [indx,contrast_all1];
contrast_ttl= zeros(numel(t),1);
trial_start = zeros(numel(t),1);
response1 = zeros(numel(t),1);
for i = 1:numel(allPD2)
    contrast_ttl([contrast_ttl1(i,1)+1:contrast_ttl1(i,1)+5],:) = repmat(contrast_ttl1(i,2)*4,5,1);
    trial_start([contrast_ttl1(i,1)+1:contrast_ttl1(i,1)+5],:) = 4;
    if not(response(i)== 0)
        delay = round(response_time(i)*35);
        response1([contrast_ttl1(i,1)+1+delay:contrast_ttl1(i,1)+5+delay],:) = repmat(response(i)*4,5,1);
    end
end
%%
task_label1 = T_all.label;
task_label = repmat("repeat",numel(indx),1);
task_label(norepeat_indx) = task_label1;
%
trial_onset2 = nan(numel(task_label),1);
trial_onset2(norepeat_indx) = T_all.wheel_onset;
trial_onset2 = trial_onset2+allPD2;
%
trial_correct = nan(numel(t),1);
trial_incorrect = nan(numel(t),1);
trial_miss = nan(numel(t),1);
trial_falarm = nan(numel(t),1);
trial_reject = nan(numel(t),1);
trial_repeat = nan(numel(t),1);
for i = 1:numel(task_label)
    current_label = task_label(i);
    switch current_label
        case 'correct'
            trial_correct([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_correct([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_correct([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
        case 'incorrect'
            trial_incorrect([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_incorrect([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_incorrect([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
        case 'miss'
            trial_miss([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_miss([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_miss([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
        case {'falarmL', 'falarmR'}
            trial_falarm([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_falarm([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_falarm([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
        case 'reject'
            trial_reject([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_reject([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_reject([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
        case 'repeat'
            trial_repeat([contrast_ttl1(i,1):contrast_ttl1(i,1)+5]) = 4;
            trial_repeat([contrast_ttl1(i,1)-2:contrast_ttl1(i,1)-1]) = 0;
            trial_repeat([contrast_ttl1(i,1)+6:contrast_ttl1(i,1)+7]) = 0;
    end
end
%%
hs14b =figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
plot(t,wf1(1,:)/2,'k');
hold on;
plot(t,wf1(2,:)/2+5,'k');
hold on;
plot(t,wf1(5,:)/2+10,'k');
hold on;

plot(t,vela/500+25,'b');
hold on;
scatter(trial_onset2,25,6,'r','filled');
hold on;
plot(t,contrast_ttl+45,'k');
hold on;
plot(t,response1+35,'k');


label_all = {'correct','incorrect','miss','reject','falarm','repeat'};
color_all = {'k','y','r','g','b','m'};
hold on;
plot(t,trial_start+55,'k');
hold on;
plot(t,trial_correct+55,'k');
hold on;
plot(t,trial_incorrect+55,'y');
hold on;
plot(t,trial_miss+55,'r');
hold on;
plot(t,trial_reject+55,'g');
hold on;
plot(t,trial_falarm+55,'b');
hold on;
plot(t,trial_repeat+55,'m');

xline(t(indx),'--k');

yticks([0,5,10,25,35,45,55]);
yticklabels({'V1-L','V1-R','SSp-ul',...
    'Wheel','Reward','Contrast','Trial-start'});
annotation('textbox', [0, 0.9, 0, 0], 'string', 'correct','color',color_all{1});
annotation('textbox', [0, 0.85, 0, 0], 'string', 'incorrect','color',color_all{2});
annotation('textbox', [0, 0.8, 0, 0], 'string', 'miss','color',color_all{3});
annotation('textbox', [0, 0.75, 0, 0], 'string', 'reject','color',color_all{4});
annotation('textbox', [0, 0.7, 0, 0], 'string', 'falarm','color',color_all{5});
annotation('textbox', [0, 0.65, 0, 0], 'string', 'repeat','color',color_all{6});
tt1 = 1920; tt2 = 2020;
xlim([tt1, tt2]);
print(hs14b, fullfile(save_folder,['FigS14b_' mn '_session' num2str(kk) '_example_'...
    num2str(tt1) 'to' num2str(tt2) 's']),'-dpdf', '-bestfit', '-painters');