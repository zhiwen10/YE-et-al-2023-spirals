function hs8l = plotWaveIndexAmp(data_folder,save_folder)
%%
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all.mat'));
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_amp_all.mat'));
vxy_all = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
angle_all = angle(vxy_all);
sync_amp_all = abs(vxy_all);
osci_amp_all = cat(3,traceAmpt_MO_left,traceAmpt_MO_right,traceAmpt_SSp_left,traceAmpt_SSp_right);
edges = [0:0.00125:0.025];
sync_amp_bins = nan(numel(edges)-1,15);
for kk = 1:15
    angle_temp = squeeze(angle_all(kk,:,4));
    sync_amp_temp = squeeze(sync_amp_all(kk,:,4));
    osci_amp_temp = squeeze(osci_amp_all(kk,:,4));
    osci_amp_temp = osci_amp_temp(1:end-1);
    %%  
    [N,edges] = histcounts(osci_amp_temp,edges);
    for i = 1:numel(edges)-1
        clear a 
        a = find(osci_amp_temp>=edges(i) & osci_amp_temp<edges(i+1));
        if not(isempty(a))       
            sync_amp_bins(i,kk) = mean(sync_amp_temp(a));
        end
    end
end
sync_amp_bins_mean = mean(sync_amp_bins,2,'omitnan');
sync_amp_bins_sem = std(sync_amp_bins,[],2,'omitnan')./sqrt(15);
hs8l = figure('Renderer', 'painters', 'Position', [50 50 250 350]);
errorbar(edges(1,1:8)*100,sync_amp_bins_mean(1:8),sync_amp_bins_sem(1:8),'k');
xticks([0,0.5,1]);
xticklabels({'0','0.5','1'});
xlabel({'2-8Hz amp(dF/F, %)'});
ylabel('Plane wave index');
xlim([0,1]);
ylim([0.2,0.6]);
yticks([0.2:0.1:0.6]);
yticklabels({'0.2','0.3','0.4','0.5','0.6'});
print(hs8l, fullfile(save_folder,'FigS8l_planewave_index_vs_amp.pdf'),...
    '-dpdf', '-bestfit', '-painters');