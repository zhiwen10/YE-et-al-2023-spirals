function [wf_all,filt_all,phase_all,amp_all, contrast_all] = getTrialTraceTask3(data_folder,T1,win,trialWin, freq)
downscale = 16;
Visp(1,:) = [810,280]; % left hemisphere V1
Visp(2,:) = [810,880]; % right hemisphere V1
Visp(3,:) = [900,280]; % left hemisphere V1 edge
Visp(4,:) = [900,880]; % right hemisphere V1 edge
Visp(5,:) = [512,850]; % right SSp-ul
Visp(6,:) = [340,705]; % right ALM
Visp_ds = round(Visp./downscale);
wf_all = [];
filt_all = [];
phase_all = [];
amp_all = [];
contrast_all = [];
%%
nStart = -trialWin(1)*35;
nEnd = trialWin(2)*35;
%%
for kk = 1:size(T1,1)
    %%
    clear contrasts Va Va1 wf wf1 pd tlTimes tt fs1 wf1a wf1b
    %%
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
    allPD2 = getPhotodiodeTime(session_root,win);
    %% load atlas registration
    load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
    sizeTemplate = [1320,1140];
    Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
    mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
    mimgt = mimgt(1:downscale:end,1:downscale:end); 
    Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
    dV1 = double(V(1:50,:));
    %% get all traces for 4 pixels
    wf1 = [];
    for j = 1:6
        wfa = squeeze(Ut1(Visp_ds(j,1),Visp_ds(j,2),:))'*dV1;
        wfa = wfa./mimgt(Visp_ds(j,1),Visp_ds(j,2));
        wf1 = [wf1;wfa];
    end
    wf1 = wf1';
    %%
    trace1_demean = wf1-mean(wf1,1);
    trace1_demean = double(trace1_demean);
    %%
    Fs = 35;
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    traceFilt = filtfilt(f1,f2,double(trace1_demean));
    traceHilbert =hilbert(traceFilt);
    tracePhase = angle(traceHilbert);
    traceAmp = abs(traceHilbert);
    %%
    for i = 1:numel(allPD2)
        indx = find(t-allPD2(i)>0,1,'first');
        wf1b(:,:,i) = wf1(indx-nStart:indx+nEnd,:);
        filt1b(:,:,i) = traceFilt(indx-nStart:indx+nEnd,:);
        phase1b(:,:,i) = tracePhase(indx-nStart:indx+nEnd,:);
        amp1b(:,:,i) = traceAmp(indx-nStart:indx+nEnd,:);
    end  
    %%
    ntrial = numel(block.events.endTrialValues);
    norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
    response = block.events.responseValues(1:ntrial);
    norepeat_indx = (norepeatValues==1);
    %%
    left_contrast = block.events.contrastLeftValues(1:ntrial);
    right_contrast = block.events.contrastRightValues(1:ntrial);
    left_contrast_unique = left_contrast(norepeat_indx);
    right_contrast_unique = right_contrast(norepeat_indx);
    contrasts = [left_contrast_unique; right_contrast_unique]';
    %%
    wf1b = wf1b(:,:,1:ntrial);
    wf1b = wf1b(:,:,norepeat_indx); 
    filt1b = filt1b(:,:,1:ntrial);
    filt1b = filt1b(:,:,norepeat_indx); 
    phase1b = phase1b(:,:,1:ntrial);
    phase1b = phase1b(:,:,norepeat_indx); 
    amp1b = amp1b(:,:,1:ntrial);
    amp1b = amp1b(:,:,norepeat_indx); 
    %% concatenate
    wf_all = cat(3,wf_all,wf1b);    
    phase_all = cat(3,phase_all,phase1b);
    filt_all = cat(3,filt_all,filt1b);
    amp_all = cat(3,amp_all,amp1b);
    contrast_all = [contrast_all;contrasts];
end