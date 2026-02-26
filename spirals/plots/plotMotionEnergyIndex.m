function hs8j = plotMotionEnergyIndex(data_folder,save_folder)
load(fullfile(data_folder,'spirals','spirals_index',...
    'energy_sync_bin_0_0.01_all2.mat'));
%%
hs8j = figure;
for i = 1:9
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

for i = 10:13
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
xlabel('Motion energy');
ylabel('sync index');
ylim([0,1]);
subplot(1,3,2)
xlabel('Motion energy');
ylabel('spirality index');
ylim([0,1]);
subplot(1,3,3)
xlabel('Motion energy');
ylabel('order index');
ylim([0,1]);
%%
print(hs8j, fullfile(save_folder,'Figs8j_motion_vs_orderness'),...
    '-dpdf', '-bestfit', '-painters');