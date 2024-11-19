function getMeanMapSession(data_folder,save_folder)
%% define params
save_folder1 = fullfile(save_folder, 'individual');
win = [0,5000];
trialWin = [-4,4];
downscale = 16;
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels = {"correct", "incorrect","miss"};
contrasts = [1,0.5,0.25,0.125,0.06,0,0,0,0,0,0;...
    0,0,0,0,0,1,0.5,0.25,0.125,0.06,0];
%%
%for each mouse and trial type, loop through all sessions, and concanetate
%all trials
for i = 1:4
    %%
    mn = fnames{i};
    T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
    T1 = T_session(T_session.label == "task",:);
    T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
    %
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome.mat']));
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_trial_ID.mat']));
    %%
    for type_i = 1:3
        %%
        label = labels{type_i};
        %
        clear wf_all wf
        correct_trials = sum(T_all.label== label);
        wf_all = zeros(83,72,141,correct_trials);
        contrast_all = [];
        count_all= 0;
        count1 = 1;
        for kk = 1:size(T1,1)
            clear Va ta wf tt indx
            %
            mn = char(T1.MouseID{kk});
            tda = T1.date(kk);  
            en = T1.folder(kk); 
            block_en = T1.block_folder(kk); 
            td = datestr(tda,'yyyy-mm-dd');
            tdb = datestr(tda,'yyyymmdd');
            % load data and block 
            fname = [mn '_' tdb '_' num2str(en)];
            session_root = fullfile(data_folder,'task','task_svd',fname);
            [U,V,t,mimg] = loadUVt2(session_root);                                 % load U,V, t
            dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
            expDir = dir(fullfile(session_root,'*_Block.mat'));
            load(fullfile(expDir.folder,expDir.name));
            % get photodiode time
            allPD2 = getPhotodiodeTime(session_root,win);
            % load atlas registration
            load(fullfile(data_folder,'task','rfmap',[fname '.mat']));
            sizeTemplate = [1320,1140];
            Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
            mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
            mimgt = mimgt(1:downscale:end,1:downscale:end); 
            Ut1 = Ut(1:downscale:end,1:downscale:end,1:50)./mimgt; 
            dV1 = double(dV(1:50,:));
            %
            ntrial = numel(block.events.endTrialValues);
            norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
            response = block.events.responseValues(1:ntrial);
            norepeat_indx = (norepeatValues==1);
            allPD2 = allPD2(1:ntrial);
            allPD2 = allPD2(norepeat_indx);
            % get time index for all events
            clear indx
            trial_current = trial_all(trial_all.session == kk ,:);
            allPD3 = allPD2(trial_current.label == label);
            if not(isempty(allPD3))
                for i = 1:numel(allPD3)
                    ta = t-allPD3(i);
                    indx(i,1) = find(ta>0, 1,'first');
                    tt(i,:) = allPD3(i)-2:1/35:allPD3(i)+2;
                    Va(:,:,i) = dV1(:,indx(i,1)-70:indx(i,1)+70);
                end
                % image sequence for mean trace map
                % only look at non-repeat trials
                Va = reshape(Va,size(Va,1),size(Va,2)*size(Va,3));
                wf = reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3))*Va;
                wf = reshape(wf, size(Ut1,1),size(Ut1,2),141,numel(allPD3));  
            else
                wf = [];
            end
            % concatenate
            count2 = count1+numel(allPD3)-1;
            wf_all(:,:,:,count1:count2) = wf;
            count1 = count1+numel(allPD3);
            count_all = count_all++numel(allPD3);
        end
        %%
        Ta = T_all(T_all.label== label,:);
        trace_correct_mean = nan([size(wf_all,[1,2,3]),11]);
        for ii = 1:11
            clear trace   
            index = (Ta.left_contrast == contrasts(1,ii) & ...
                Ta.right_contrast == contrasts(2,ii));
            trace = squeeze(mean(wf_all(:,:,:,index),4));
            trace_correct_mean(:,:,:,ii) = trace;
        end
        %%
        save(fullfile(save_folder1,[mn '_mean_map_' char(label) '.mat']),'contrasts','trace_correct_mean');
    end
end