function h5gh = plotHitRateSlow(data_folder,save_folder)
%%
fname = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
hit_rate_high = nan(5,4);
miss_rate_high = nan(5,4);
hit_rate_low = nan(5,4);
miss_rate_low = nan(5,4);
rt_high = nan(5,4);
rt_low = nan(5,4);
freq = [0.05 2];
%%
for m = 1:4
    clear T_session T_all wf_all wft T_high T_low
    mn = fname{m};
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_freq_to' num2str(freq(2)) 'Hz']));
    wf_all = wf_all*100;
    %%
    amp = mean(squeeze(amp_all(1:141,1,:)),1)';
    p1 = prctile(amp,25);
    p2 = prctile(amp,75);
    T_high = T_all(amp<p1,:);
    T_low = T_all(amp>p2,:);
    %% right hemisphere   
    [T_high_ratio,rt_high_median] = sort_ratio_by_contrast2(T_high);
    [T_low_ratio,rt_low_median] = sort_ratio_by_contrast2(T_low);
    %%
    hit_rate_high(:,m) = T_high_ratio{:,3};
    miss_rate_high(:,m) = T_high_ratio{:,2};
    rt_high(:,m) = rt_high_median;
    hit_rate_low(:,m) = T_low_ratio{:,3};
    miss_rate_low(:,m) = T_low_ratio{:,2};
    rt_low(:,m) = rt_low_median;
end
%%
subj = [1,2,3,4];
subjs = repmat(subj,[5,1]);
subject = [subjs;subjs];
hit_rates = [hit_rate_high;hit_rate_low];
miss_rates = [miss_rate_high; miss_rate_low];
rt = [rt_high;rt_low];
contrast_matrix1 = [6.25;12.5;25;50;100];
contrast_matrix2 = repmat(contrast_matrix1,[1,4]);
contrasts = [contrast_matrix2;contrast_matrix2];
voltage_label1 = ones(5,4);
voltage_label2 = ones(5,4)+1;
voltage_label = [voltage_label1;voltage_label2];
%%
subject = categorical(subject(:));
rt = rt(:);
hit_rates = hit_rates(:);    
miss_rates = miss_rates(:);
contrasts = categorical(contrasts(:));
voltage_label2 = categorical(voltage_label(:));
voltage_label3 = voltage_label(:);
label = strings(40,1);
label(voltage_label3 ==1) = 'high';
label(voltage_label3 ==2) = 'low';
%%
pp_hit = anovan(hit_rates,{contrasts, voltage_label2,subject},...
    'model',2,'random',3,'varnames',{'contrast','label','subj'});
%%
hit_rate_change = hit_rate_high-hit_rate_low;
miss_rate_change = miss_rate_high-miss_rate_low;
hit_rate_change_mean = mean(hit_rate_change,2);
hit_rate_change_sem = std(hit_rate_change,[],2)./sqrt(6);
% % one way repeated anova (wrong test, one group is significantly differently than all the rest)
% tbl = array2table(hit_rate_change','VariableNames',{'A','B','C','D','E'});
% withinDesign = table([1 2 3 4 5]', 'VariableNames', {'Contrast'});
% withinDesign.Contrast = categorical(withinDesign.Contrast);
% % Fit the linear model
% rm = fitrm(tbl,'A-E~1','WithinDesign', withinDesign);
% ranova(rm, 'WithinModel', 'Contrast')
%%
h5gh =figure('Renderer', 'painters', 'Position', [100 100 600 300]);
subplot(1,2,1);
subjects = ones(4,1);
for i = 1:5
    scatter(subjects*4*i,hit_rate_high(i,:),8,'b');
    hold on;
    scatter(subjects*4*i-2,hit_rate_low(i,:),8,'m');
    hold on;
    plot([subjects*4*i,subjects*4*i-2]',[hit_rate_high(i,:);hit_rate_low(i,:)],'k');
end
xticks([3,7,11,15,19]);
xticklabels({'6%','12.5%','25%','50%','100%'});
xlabel('Contrast');
ylabel('Hit rate');
ylim([0,1.0]);

subplot(1,2,2);
for i = 1:5
    scatter(subjects*4*i-1,hit_rate_change(i,:),8,'k');
    hold on;
end
errorbar(4*[1:5]-1, hit_rate_change_mean,hit_rate_change_sem,'color','r','linewidth',2);
yline(0,'--k');
ylim([-0.3,0.3]);
xticks([3,7,11,15,19]);
xticklabels({'6%','12.5%','25%','50%','100%'});
xlabel('Contrast');
ylabel('Hit rate change');
%%
print(h5gh, fullfile(save_folder,'Fig5gh_hit_ratio_005_2Hz'),'-dpdf', '-bestfit', '-painters');
