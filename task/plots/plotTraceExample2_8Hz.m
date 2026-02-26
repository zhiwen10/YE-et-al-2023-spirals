function h5i = plotTraceExample2_8Hz(data_folder,save_folder)
%%
mn = 'ZYE_0091';
trialWin = [-4,4];
freq = [2,8];
%%
load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
load(fullfile(data_folder,'task','task_outcome',[mn '_task_freq_to' num2str(freq(2)) 'Hz']));
wf_all = wf_all*100;
%%
t1 = trialWin(1):1/35:trialWin(2);
stimOn = (numel(t1)-1)/2+1;
amp1 = squeeze(amp_all(:,1,:));
amp2 = sum(amp1(123:141,:),1);
p1 = prctile(amp2,90);
indx = (amp2>p1);
contrast_all1 = contrast_all(indx,:);
T_all1 = T_all(indx,:);
wf_all1 = squeeze(wf_all(:,1,indx));
%%
trace1_demean = wf_all1-mean(wf_all1,1);
trace1_demean = double(trace1_demean);
Fs = 35;
freq = [2,8];
[f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
traceFilt = filtfilt(f1,f2,double(trace1_demean));
traceHilbert =hilbert(traceFilt);
tracePhase = angle(traceHilbert);
traceAmp = abs(traceHilbert);
%%
index = (abs(contrast_all1(:,1)-contrast_all1(:,2)) >0 & abs(contrast_all1(:,1)-contrast_all1(:,2)) <=0.125);
traceRaw = squeeze(wf_all1(:,index));
taskR = T_all1(index,:);
label_all = taskR.label;
scale = 1;
label = {'correct','miss'};
color1 = {'k','r'};
h5i = figure('Renderer', 'painters', 'Position', [50  50 450 450]);
for m = 1:2
    clear index1 traceRaw1
    index1 = strcmp(label_all,label{m});
    traceRaw1 = traceRaw(:,index1);
    N = size(traceRaw1,2);
    traceRaw2 = traceRaw1(:,1:20);
    traceRaw1_mean = squeeze(mean(traceRaw1,2));
    
    subplot(3,2,m);
    plot(t1,squeeze(traceRaw1_mean),color1{m});
    hold on;
    xline(t1(stimOn),'--k');
    ylim([-4*scale,4*scale]);
    
    % let's only use 15-32 trials as examples, since they look good
    if m ==1
        traceRaw2 = traceRaw2(:,[2,5,6,7,9,11,12,14,19,20]);
    else
        traceRaw2 = traceRaw2(:,[1,4,5,8,11,13,14,18,19,20]);
    end
    subplot(3,2,[2+m,4+m]);
    for i = 1:10
        plot(t1,squeeze(traceRaw2(:,i))+10*i*scale,color1{m});
        hold on;
    end
    hold on;
    plot([-3.5, -3.5],[0,8],'r');
    xline(t1(stimOn),'--k');
    ylim([-10,110]);
end
print(h5i, fullfile(save_folder,'Fig5i_example_2_8Hz'),'-dpdf', '-bestfit', '-painters');