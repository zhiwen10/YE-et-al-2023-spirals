%%
h1 = figure('Renderer', 'painters', 'Position', [100 100 400 800]);
% plot instantenous frequency distribution -raw data
load('instantenous_frequency_across_sessions_V.mat');
load('instantenous_frequency_across_sessions_V_all.mat');
edges = 0:0.5:15;
for n = 1:15
    clear counts_temp1 counts_temp2
    frt_frame_spiral = frt_all{n};
    counts_temp1 = histcounts(frt_frame_spiral(:,1),edges);
    N_frame_spiral(:,n) = counts_temp1;
    ratio_frame_spiral(:,n) = N_frame_spiral(:,n)./size(frt_frame_spiral,1);
    %%
    frt_frame_all = frt_all1{n};
    counts_temp2 = histcounts(frt_frame_all(:,1),edges);
    N_frame_all(:,n) = counts_temp2;
    ratio_frame_all(:,n) = N_frame_all(:,n)./size(frt_frame_all,1);
end
ratio_frame_sprial_mean = mean(ratio_frame_spiral,2);
ratio_frame_spiral_sem = std(ratio_frame_spiral,[],2)./sqrt(15);
ratio_frame_all_mean = mean(ratio_frame_all,2);
ratio_frame_all_sem = std(ratio_frame_all,[],2)./sqrt(15);
subplot(3,2,1);
s2 = errorbar(edges(1:end-1),ratio_frame_sprial_mean,ratio_frame_spiral_sem,'r');
hold on;
s2 = errorbar(edges(1:end-1),ratio_frame_all_mean,ratio_frame_all_sem,'k');
yticks([0:0.2:0.8]);
ylim([0,0.8]);
yticklabels({[0:0.2:0.8]});
xlabel('Instantenous frequency (Hz)');
ylabel('Proportion of spirals');
set(gca, 'box', 'off');
subplot(3,2,2);
for n = 1:15
    clear frt
    frt = frt_all{n};
    for mm = 1:14
        clear frt_temp
        frt_temp = frt(frt(:,2)==mm,1);
        frt_mean(n,mm) = mean(frt_temp,"omitnan");
        frt_sem(n,mm) = std(frt_temp,[],1,"omitnan")./sqrt(length(frt_temp));
    end
end
frt_mean2 = mean(frt_mean,1,"omitnan");
count1 = sum(not(isnan(frt_mean)),1);
frt_std2 = std(frt_mean,[],1,"omitnan");
for m = 1:14
    frt_sem2(m,1) = frt_std2(m)./sqrt(count1(m));
end
errorbar([1:10]/35,frt_mean2(1:10),frt_sem2(1:10),'r');
xlim([0,0.3]);
ylim([0,2]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');

% plot instantenous frequency distribution -dv data
load('instantenous_frequency_across_sessions_dV.mat');
load('instantenous_frequency_across_sessions_dV_all.mat');
edges = 0:0.5:15;
for n = 1:15
    clear counts_temp1 counts_temp2
    frt_frame_spiral = frt_all{n};
    counts_temp1 = histcounts(frt_frame_spiral(:,1),edges);
    N_frame_spiral(:,n) = counts_temp1;
    ratio_frame_spiral(:,n) = N_frame_spiral(:,n)./size(frt_frame_spiral,1);
    %%
    frt_frame_all = frt_all1{n};
    counts_temp2 = histcounts(frt_frame_all(:,1),edges);
    N_frame_all(:,n) = counts_temp2;
    ratio_frame_all(:,n) = N_frame_all(:,n)./size(frt_frame_all,1);
end
ratio_frame_sprial_mean = mean(ratio_frame_spiral,2);
ratio_frame_spiral_sem = std(ratio_frame_spiral,[],2)./sqrt(15);
ratio_frame_all_mean = mean(ratio_frame_all,2);
ratio_frame_all_sem = std(ratio_frame_all,[],2)./sqrt(15);
subplot(3,2,3);
s2 = errorbar(edges(1:end-1),ratio_frame_sprial_mean,ratio_frame_spiral_sem,'r');
hold on;
s2 = errorbar(edges(1:end-1),ratio_frame_all_mean,ratio_frame_all_sem,'k');
[ma,mindex] = max(ratio_frame_sprial_mean);
s3 = xline(edges(mindex),'--');
text(edges(mindex),0.08,['maxFreq = ' num2str(edges(mindex)) 'Hz']);
yticks([0:0.02:0.1]);
ylim([0,0.1]);
yticklabels({[0:0.02:0.1]});
xlabel('Instantenous frequency (Hz)');
ylabel('Proportion of spirals');
set(gca, 'box', 'off');
subplot(3,2,4);
for n = 1:15
    clear frt
    frt = frt_all{n};
    for mm = 1:14
        clear frt_temp
        frt_temp = frt(frt(:,2)==mm,1);
        frt_mean(n,mm) = mean(frt_temp,"omitnan");
        frt_sem(n,mm) = std(frt_temp,[],1,"omitnan")./sqrt(length(frt_temp));
    end
end
frt_mean2 = mean(frt_mean,1,"omitnan");
count1 = sum(not(isnan(frt_mean)),1);
frt_std2 = std(frt_mean,[],1,"omitnan");
for m = 1:14
    frt_sem2(m,1) = frt_std2(m)./sqrt(count1(m));
end
errorbar([1:10]/35,frt_mean2(1:10),frt_sem2(1:10),'r');
xlim([0,0.3]);
ylim([4,8]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');

% plot instantenous frequency distribution -filt data
load('instantenous_frequency_across_sessions_dV_filt.mat');
load('instantenous_frequency_across_sessions_dV_filt_all.mat');
edges = 0:0.5:15;
for n = 1:15
    clear counts_temp1 counts_temp2
    frt_frame_spiral = frt_all{n};
    counts_temp1 = histcounts(frt_frame_spiral(:,1),edges);
    N_frame_spiral(:,n) = counts_temp1;
    ratio_frame_spiral(:,n) = N_frame_spiral(:,n)./size(frt_frame_spiral,1);
    %%
    frt_frame_all = frt_all1{n};
    counts_temp2 = histcounts(frt_frame_all(:,1),edges);
    N_frame_all(:,n) = counts_temp2;
    ratio_frame_all(:,n) = N_frame_all(:,n)./size(frt_frame_all,1);
end
ratio_frame_sprial_mean = mean(ratio_frame_spiral,2);
ratio_frame_spiral_sem = std(ratio_frame_spiral,[],2)./sqrt(15);
ratio_frame_all_mean = mean(ratio_frame_all,2);
ratio_frame_all_sem = std(ratio_frame_all,[],2)./sqrt(15);
subplot(3,2,5);
s2 = errorbar(edges(1:end-1),ratio_frame_sprial_mean,ratio_frame_spiral_sem,'r');
hold on;
s2 = errorbar(edges(1:end-1),ratio_frame_all_mean,ratio_frame_all_sem,'k');
[ma,mindex] = max(ratio_frame_sprial_mean);
s3 = xline(edges(mindex),'--');
text(edges(mindex),0.08,['maxFreq = ' num2str(edges(mindex)) 'Hz']);
yticks([0:0.1:0.3]);
ylim([0,0.3]);
yticklabels({[0:0.1:0.3]});
xlabel('Instantenous frequency (Hz)');
ylabel('Proportion of spirals');
set(gca, 'box', 'off');
subplot(3,2,6);
for n = 1:15
    clear frt
    frt = frt_all{n};
    for mm = 1:14
        clear frt_temp
        frt_temp = frt(frt(:,2)==mm,1);
        frt_mean(n,mm) = mean(frt_temp,"omitnan");
        frt_sem(n,mm) = std(frt_temp,[],1,"omitnan")./sqrt(length(frt_temp));
    end
end
frt_mean2 = mean(frt_mean,1,"omitnan");
count1 = sum(not(isnan(frt_mean)),1);
frt_std2 = std(frt_mean,[],1,"omitnan");
for m = 1:14
    frt_sem2(m,1) = frt_std2(m)./sqrt(count1(m));
end
errorbar([1:10]/35,frt_mean2(1:10),frt_sem2(1:10),'r');
xlim([0,0.3]);
ylim([4,5]);
set(gca, 'box', 'off');
xlabel('Spiral duration (s)');
ylabel('Instantenous frequency (Hz)');
print(h1, 'Instantenous frequency_session_summary.pdf','-dpdf', '-bestfit', '-painters');