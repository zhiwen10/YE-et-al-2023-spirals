function getMeanMapsAll(data_folder,save_folder)
%%
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
labels ={'correct','incorrect', 'miss'};
freq = [2,8];
Fs = 35;
%%
trace_mean_all = zeros(83,72,141,3);
traceFilt_all = zeros(83,72,141,3);
tracePhase_all = zeros(83,72,141,3);
for kk = 1:3
    label = labels{kk};
    trace_correct_mean_all = [];
    for ii = 1:4
        mn = fnames{ii};
        load(fullfile(data_folder,'task','task_mean_maps','individual',[mn '_mean_map_' label '.mat']));
        trace_correct_mean_all = cat(5,trace_correct_mean_all,trace_correct_mean);
    end
    %%
    trace_correct_mean_all1 = squeeze(mean(trace_correct_mean_all,5,'omitnan'));
    cts = [1:3];
    trace_mean = squeeze(mean(trace_correct_mean_all1(:,:,:,cts),4,'omitnan')); % only for first 4 contrasts
    trace_mean2 = reshape(trace_mean,[size(trace_mean,1)*size(trace_mean,2),size(trace_mean,3)]);
    %% take care of nan pixels after registration before filting 
    indx1 = not(isnan(trace_mean2(:,1)));
    trace_mean4 = trace_mean2(indx1,:);
    meanTrace = trace_mean4 -mean(trace_mean4,2);
    meanTrace = meanTrace';  
    [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
    traceFilt = filtfilt(f1,f2,meanTrace);
    traceHilbert =hilbert(traceFilt);
    tracePhase = angle(traceHilbert);
    tracePhase2 = nan(size(trace_mean2));
    traceFilt2 = nan(size(trace_mean2));
    tracePhase2(indx1,:) = tracePhase';
    traceFilt2(indx1,:) = traceFilt';
    %%
    trace_mean2 = trace_mean2-mean(trace_mean2,2);
    trace_mean3 = reshape(trace_mean2,size(trace_mean,1),size(trace_mean,2),size(trace_mean2,2));
    traceFilt3 = reshape(traceFilt2,size(trace_mean,1),size(trace_mean,2),size(traceFilt2,2));
    tracePhase3 = reshape(tracePhase2,size(trace_mean,1),size(trace_mean,2),size(tracePhase2,2));
    %%
    trace_mean_all(:,:,:,kk) = trace_mean3;
    traceFilt_all(:,:,:,kk) = traceFilt3;
    tracePhase_all(:,:,:,kk) = tracePhase3;
end
%%
save(fullfile(save_folder,'task_mean_maps_all_mice.mat'),'trace_mean_all','traceFilt_all','tracePhase_all');