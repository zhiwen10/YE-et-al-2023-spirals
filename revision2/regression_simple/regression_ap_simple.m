githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));       % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% Figure 3c
save_folder = fullfile(data_folder, 'spirals_mirror','regression_ap_simple');
getReducedRankRegressionAP_simple(T,data_folder,save_folder);                     % save regression coeffs, kernels and R2 from AP regression
%%
save_folder = fullfile(data_folder, 'spirals_mirror','regression_hemi_simple');
getReducedRankRegressionHEMI_simple(T,data_folder,save_folder);                     % save regression coeffs, kernels and R2 from AP regression
%%
save_folder = fullfile(data_folder, 'spirals_mirror','regression_hemi_simple2');
getReducedRankRegressionHEMI_simple_invert(T,data_folder,save_folder);                     % save regression coeffs, kernels and R2 from AP regression
%%
save_folder1 = fullfile(data_folder, 'spirals_mirror');
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    % load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en) '-AP_simple.mat'];
    load(fullfile(save_folder1,'regression_ap_simple',fname));
    var_all_ap(kk,:) =  explained_var;    
    %%
    fname2 = [mn '_' tdb '_' num2str(en) '-HEMI_simple.mat'];
    load(fullfile(save_folder1,'regression_hemi_simple',fname2));
    var_all_hemi(kk,:) =  explained_var;
end
%%
h3c = figure('Renderer', 'painters', 'Position', [100 100 700 300]);
subplot(1,2,1);
mean_var_ap_all = mean(var_all_ap,1);
sem_var_ap_all = std(var_all_ap,[],1)./sqrt(15);
shadedErrorBar([-10:2:10]/35, mean_var_ap_all, sem_var_ap_all, 'lineprops', '-r');
hold on;
xline(0,'k--');
subplot(1,2,2);
mean_var_hemi_all = mean(var_all_hemi,1);
sem_var_hemi_all = std(var_all_hemi,[],1)./sqrt(15);
shadedErrorBar([-10:2:10]/35, mean_var_hemi_all, sem_var_hemi_all, 'lineprops', '-g');
hold on;
xline(0,'k--');
%%
[h1,p1] = ttest(var_all_hemi(:,4),var_all_hemi(:,8));