function h5f = plotMeanTraceExampleSlow(data_folder,save_folder)
%%
h5f = figure('Renderer', 'painters', 'Position', [50 50 950 250]);
mn = 'ZYE_0091';
T_session = readtable(fullfile(data_folder,'task','sessions',[mn '.xlsx']));
T1 = T_session(T_session.label == "task",:);
T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
%
trialWin = [-4,4];
load(fullfile(data_folder,'task','trial_trace',[mn '_task']));
wf_all = wf_all*100;
t1 = trialWin(1):1/35:trialWin(2);
stimOn = (numel(t1)-1)/2+1;
contrasts = [0.06,0.125,0.25,0.5,1;0,0,0,0,0];
color1 = {'k','r'};
label = {'correct','incorrect','miss'};
for ii = 1:5
    clear trace1 traceRaw taskR
    index = (contrast_all(:,1) == contrasts(1,ii) & contrast_all(:,2) == contrasts(2,ii)); % 1 is left, 2 is right
    trace1 = squeeze(wf_all(:,:,index));
    traceRaw = permute(trace1,[2,1,3]);
    taskR = T_all(index,:);
    scale = 1;
    label_all = taskR.label;
    for m = 1:2
        subplot(1,10,(ii-1)*2+m);
        clear index1 traceRaw1
        if m==1
            label = 'correct';
        else
            label = {'miss'};
        end
        index1 = ismember(label_all,label);
        traceRaw1 = traceRaw(:,:,index1);
        traceRaw1_mean = mean(squeeze(traceRaw1(:,1,:)),2);

        plot(t1,squeeze(traceRaw1_mean),color1{m});
        hold on;
        xline(t1(stimOn),'--k');
        ylim([-4*scale,4*scale]);
        if m==1
            title(['C' num2str(floor(contrasts(1,ii)*100))]);
        end
    end
end
%%
print(h5f, fullfile(save_folder,'Fig5f_mean_trace_example_slow'), ...
    '-dpdf', '-bestfit', '-painters');