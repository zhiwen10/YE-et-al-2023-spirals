function hs14cd = plotMeanTraceAcrossContrasts(data_folder,save_folder)
%%
fname = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
freq = [2 8];
cs = [0.06,0.125,0.25,0.5,1];
t1 = -0.5:1/35:0.5;
pixel_index = 124:159;
for m = 1:4
    clear T_session T_all wf_all wft T_high T_low index phase1
    mn = fname{m};
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_outcome']));
    load(fullfile(data_folder,'task','task_outcome',[mn '_task_freq_to' num2str(freq(2)) 'Hz']));
    wf_all = wf_all*100;
    %%
    for kk = 1:5
        clear index1 index2 trace_all1 trace_all2 trace_all3
        index1 = (T_all.label == "correct" & (T_all.left_contrast-T_all.right_contrast)== cs(kk));
        trace_all1_conrta = squeeze(wf_all(:,2,index1));
        index2 = (T_all.label == "correct" & (T_all.right_contrast-T_all.left_contrast)== cs(kk));
        trace_all2_conrta = squeeze(wf_all(:,1,index2));
        trace_all3_conrta = [trace_all1_conrta,trace_all2_conrta];
        trace_task_contra(:,kk,m) = mean(trace_all3_conrta,2);
        %%
        trace_all1_ipsi = squeeze(wf_all(:,1,index1));
        trace_all2_ipsi = squeeze(wf_all(:,2,index2));
        trace_all3_ipsi = [trace_all1_ipsi,trace_all2_ipsi];
        trace_task_ipsi(:,kk,m) = mean(trace_all3_ipsi,2);
    end
end
%%
for m = 1:4
    clear T_session T_all wf_all wft T_high T_low
    mn = fname{m};
    load(fullfile(data_folder,'task','task_outcome',[mn '_passive_freq_to' num2str(freq(2)) 'Hz']));
    wf_all = wf_all*100;
    %%
    for kk = 1:5
        clear index trace_all1
        index = ((contrast_all(:,1)-contrast_all(:,2))== cs(kk));
        trace_all_contra = squeeze(wf_all(:,2,index));
        trace_passive_contra(:,kk,m) = mean(trace_all_contra,2);
        %%
        trace_all_ipsi = squeeze(wf_all(:,1,index));
        trace_passive_ipsi(:,kk,m) = mean(trace_all_ipsi,2);        
    end
end
%%
color1 = cbrewer2('seq','Reds',5);
color2 = cbrewer2('seq','Greys',5);
hs14cd =figure('Renderer', 'painters', 'Position', [100 100 900 400]);
subplot(2,4,1);
for i = 1:5
    trace_temp = squeeze(trace_passive_contra(:,i,:));
    trace_passive_contra_mean1 = mean(trace_temp,2);
    trace_passive_contra_sem = std(trace_temp,[],2)./sqrt(4);
    passive_contra_mean(i,1) = trace_passive_contra_mean1(146)-trace_passive_contra_mean1(141);
    passive_contra_sem(i,1) = trace_passive_contra_sem(146);
    plot(t1,trace_passive_contra_mean1(pixel_index),'color',color1(i,:),'LineWidth',2);
    hold on;
end
xline(0.14);
ylim([-2,4]);
xlabel('Time (s)');
ylabel('dF/F (%)');

subplot(2,4,2);
for i = 1:5
    trace_temp = squeeze(trace_passive_ipsi(:,i,:));
    trace_passive_ipsi_mean1 = mean(trace_temp,2);
    trace_passive_ipsi_sem = std(trace_temp,[],2)./sqrt(4);
    passive_ipsi_mean(i,1) = trace_passive_ipsi_mean1(146)-trace_passive_ipsi_mean1(141);
    passive_ipsi_sem(i,1) = trace_passive_ipsi_sem(146);
    plot(t1,trace_passive_ipsi_mean1(pixel_index),'color',color2(i,:),'LineWidth',2);
    hold on;
end
xline(0.14);
ylim([-2,4]);
xlabel('Time (s)');
ylabel('dF/F (%)');

subplot(2,4,3);
errorbar(1:5,passive_contra_mean,passive_contra_sem,'r');
hold on;
errorbar(1:5,passive_ipsi_mean,passive_ipsi_sem,'k');
xlim([1,5]);
ylim([-1,3]);
xticks([1:5]);
xticklabels({'6','12','25','50','100'});
xlabel('Contrasts(%)');
ylabel('dF/F (%)');

ax4 = subplot(2,4,4);
imagesc([1:5]);
axis off;
colormap(ax4,color2);
cb = colorbar;
cb.Ticks = 1:5;
cb.TickLabels = num2cell([6,12,25,50,100]);

subplot(2,4,5);
for i = 1:5
    trace_temp = squeeze(trace_task_contra(:,i,:));
    trace_task_conrta_mean1 = mean(trace_temp,2);
    trace_task_contra_sem = std(trace_temp,[],2)./sqrt(4);
    task_contra_mean(i,1) = trace_task_conrta_mean1(146)-trace_task_conrta_mean1(141);
    task_contra_sem(i,1) = trace_task_contra_sem(146);
    plot(t1,trace_task_conrta_mean1(pixel_index),'color',color1(i,:),'LineWidth',2);
    hold on;
end
xline(0.14);
ylim([-2,4]);
xlabel('Time (s)');
ylabel('dF/F (%)');
subplot(2,4,6);
for i = 1:5
    trace_temp = squeeze(trace_task_ipsi(:,i,:));
    trace_task_ipsi_mean1 = mean(trace_temp,2);
    trace_task_ipsi_sem = std(trace_temp,[],2)./sqrt(4);
    task_ipsi_mean(i,1) = trace_task_ipsi_mean1(146)-trace_task_ipsi_mean1(141);
    task_ipsi_sem(i,1) = trace_task_ipsi_sem(146);
    plot(t1,trace_task_ipsi_mean1(pixel_index),'color',color2(i,:),'LineWidth',2);
    hold on;
end
xline(0.14);
ylim([-2,4]);
xlabel('Time (s)');
ylabel('dF/F (%)');

subplot(2,4,7);
% task_mean1 = [flipud(task_ipsi_mean);task_contra_mean];
% task_sem1 = [flipud(task_ipsi_sem);task_contra_sem];
errorbar(1:5,task_contra_mean,task_contra_sem,'r');
hold on;
errorbar(1:5,task_ipsi_mean,task_ipsi_sem,'k');
xlim([1,5]);
ylim([-1,3]);
xticks([1:5]);
xticklabels({'6','12','25','50','100'});
xlabel('Contrasts(%)');
ylabel('dF/F (%)');

ax8 = subplot(2,4,8);
imagesc([1:5]);
axis off;
colormap(ax8,color1);
cb1 = colorbar;
cb1.Ticks = 1:5;
cb1.TickLabels = num2cell([6,12,25,50,100]);
print(hs14cd, fullfile(save_folder,'FigS14cd_mean_trace'),'-dpdf', '-bestfit', '-painters');

