function h3c = plotVarianceExplained(T,data_folder,save_folder)
%%
ffolder1 = fullfile(data_folder,'spirals_mirror','regression_ap');
ffolder2 = fullfile(data_folder,'spirals_mirror','regression_hemi');
%%
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    % load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en) '-AP.mat'];
    load(fullfile(ffolder1,fname),'explained_var5');
    var_all_ap(:,kk) =  explained_var5;
end
mean_var_ap_all = mean(var_all_ap,2);
std_var_ap_all = std(var_all_ap,0,2);
sem_var_ap_all = std_var_ap_all/size(var_all_ap, 2);

for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    % load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en) '-hemi.mat'];
    load(fullfile(ffolder2,fname),'explained_var5');
    var_all_hemi(:,kk) =  explained_var5;
end
mean_var_hemi_all = mean(var_all_hemi,2);
std_var_hemi_all = std(var_all_hemi,0,2);
sem_var_hemi_all = std_var_hemi_all/size(var_all_hemi, 2);
%%
h3c = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,2,1);
shadedErrorBar(1:50, mean_var_ap_all, sem_var_ap_all, 'lineprops', '-r')
ylim([0.7,1]);
xline(16,'--k');
xlabel('Components');
ylabel('Variance Explained');
title('AP prediction');
subplot(1,2,2);
shadedErrorBar(1:50, mean_var_hemi_all, sem_var_hemi_all, 'lineprops', '-g')
ylim([0.7,1])
xline(16,'--k');
xlabel('Components');
ylabel('Variance Explained');
title('hemi prediction');
%%
print(h3c, fullfile(save_folder,'Fig3c_prediction_accuracy_rank.pdf'),...
    '-dpdf', '-bestfit', '-painters');