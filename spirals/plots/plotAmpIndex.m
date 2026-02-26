function hs8i = plotAmpIndex(data_folder,save_folder)
load(fullfile(data_folder,'spirals','spirals_index',...
    'amp_sync_bin_0_0.01_all.mat'));
%%
hs8i = figure;
for i = 1:11
    subplot(1,3,1)
    plot(edges(1:end-1),sync_sort(:,i),'k');
    hold on;
    subplot(1,3,2)
    plot(edges(1:end-1),spirality_sort(:,i),'k');
    hold on;
    subplot(1,3,3)
    plot(edges(1:end-1),index_sort(:,i),'k');
    hold on;
end

for i = 12:15
    subplot(1,3,1)
    plot(edges(1:end-1),sync_sort(:,i),'r');
    hold on;
    subplot(1,3,2)
    plot(edges(1:end-1),spirality_sort(:,i),'r');
    hold on;
    subplot(1,3,3)
    plot(edges(1:end-1),index_sort(:,i),'r');
    hold on;
end

subplot(1,3,1)
xlabel('3-6Hz amp');
ylabel('sync index');
ylim([0,1]);
subplot(1,3,2)
xlabel('3-6Hz amp');
ylabel('spirality index');
ylim([0,1]);
subplot(1,3,3)
xlabel('3-6Hz amp');
ylabel('order index');
ylim([0,1]);
%%
print(hs8i, fullfile(save_folder,'Figs8i_amp_orderness'),...
    '-dpdf', '-bestfit', '-painters');