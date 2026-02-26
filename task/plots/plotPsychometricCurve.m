function hs14a = plotPsychometricCurve(data_folder,save_folder)
load(fullfile(data_folder,'task','psychometric_curve','Task_performance_across_sessions.mat'));
%%
total_count_all = nan(6,11,4);
total_all = nan(11,4);
for i = 1:11
    for j = 1:4
        trial_count_temp = T_trial_counts_all{i,j};
        if not(isempty(trial_count_temp))
            total_count_temp = trial_count_temp{:,8}; % total trial N
            for k = 1:5 % sum left or right same contrasts
                total_count_temp2(k,1) = total_count_temp(k,1)+total_count_temp(11-k+1,1);
            end
            total_count_temp2(6,1) = total_count_temp(6,1);
            total_all_temp = sum(total_count_temp2);
            %%
            total_count_all(:,i,j) = total_count_temp2;
            total_all(i,j) = total_all_temp;
        end
    end
end
total_all2 = sum(total_all,1,'omitnan');
total_all_mean = round(mean(total_all2));
total_all_sem = round(std(total_all2)./sqrt(4));
sessionsa = not(isnan(total_all));
sessions = sum(sessionsa,1,'omitnan');
%%
total_count_mean = squeeze(sum(total_count_all,2,'omitnan'));
total_count_mean2 = mean(total_count_mean,2);
total_count_sem2 = std(total_count_mean,[],2)./sqrt(4);
total_count_mean3 = total_count_mean2/total_count_mean2(1);
total_count_mean3 = round(total_count_mean3*10)/10;
contrast1 = [100,50,25,12.5,6.25,0]';
trial_ratio_contrast = table(contrast1,total_count_mean3);
%%
fname = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
contrast = [-1,-0.5,-0.25,-0.125,-0.06,0,0.06,0.125,0.25,0.5,1]'; 
%%
hs14a =figure('Renderer', 'painters', 'Position', [100 100 800 900]);
for m = 1:4
    subplot(4,4,1+(m-1)*4);
    errorbar(contrast,right_choice_all_mean(:,m),right_choice_all_sem(:,m),'or-','MarkerEdgeColor','k','MarkerSize',4);
    ylim([0,1]);
    xticks([-1, -0.5, 0, 0.5, 1]);
    xticklabels({'-100', '-50', '0', '50', '100'});
    xlabel('Contrast (%)');
    ylabel('Left choice  probability');
    set(gca,'box','off')
    
    
    subplot(4,4,2+(m-1)*4);
    errorbar(contrast,no_go_all_mean(:,m),no_go_all_sem(:,m),'or-','MarkerEdgeColor','k','MarkerSize',4);
    ylim([0,1]);
    xticks([-1, -0.5, 0, 0.5, 1]);
    xticklabels({'-100', '-50', '0', '50', '100'});
    xlabel('Contrast (%)');
    ylabel('Miss/NoGo  probability');
    title([fname{m} ' (Trials: ' num2str(total_all2(m)) ...
        '; sessions: ' num2str(sessions(m)) ')'],'interpreter','none');
    set(gca,'box','off')
    
    subplot(4,4,3+(m-1)*4);
    errorbar(contrast,left_choice_all_mean(:,m),left_choice_all_sem(:,m),'or-','MarkerEdgeColor','k','MarkerSize',4);
    ylim([0,1]);
    xticks([-1, -0.5, 0, 0.5, 1]);
    xticklabels({'-100', '-50', '0', '50', '100'});
    xlabel('Contrast (%)');
    ylabel('Right choice probability');
    set(gca,'box','off')
    
    subplot(4,4,4+(m-1)*4);
    errorbar(contrast,rt_median_all_mean(:,m),rt_median_all_sem(:,m),'or-','MarkerEdgeColor','k','MarkerSize',4);
    ylim([0,2]);
    xticks([-1, -0.5, 0, 0.5, 1]);
    xticklabels({'-100', '-50', '0', '50', '100'});
    xlabel('Contrast (%)');
    ylabel('Reaction time (s)');
    set(gca,'box','off')
end
%%
left_choice_mean_100 = mean(left_choice_all_mean(11,:));
left_choice_sem_100 = std(left_choice_all_mean(11,:))./sqrt(4);
left_choice_mean_6 = mean(left_choice_all_mean(7,:));
left_choice_sem_6 = std(left_choice_all_mean(7,:))./sqrt(4);
[h1,p1] = ttest(left_choice_all_mean(11,:),left_choice_all_mean(7,:));

rt_mean_100 = mean(rt_median_all_mean(11,:));
rt_sem_100 = std(rt_median_all_mean(11,:))./sqrt(4);
rt_mean_6 = mean(rt_median_all_mean(7,:));
rt_sem_6 = std(rt_median_all_mean(7,:))./sqrt(4);
[h2,p2] = ttest(rt_median_all_mean(11,:),rt_median_all_mean(7,:));
%%
print(hs14a, fullfile(save_folder,['FigS14a_Psychometric_curve_all_mice.pdf']),'-dpdf', '-bestfit', '-painters');