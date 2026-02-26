%%
data_folder = 'E:\spiral_data_share\data\revision2\V_prediction_MO';
for kk = 1:15
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    load(fullfile(data_folder,[fname '_v_predict_MO.mat']),'explained_var_all');
    explained_var_mean(:,:,kk) = mean(explained_var_all,3);
end
%%
figure;
for i = 1:15
    subplot(3,5,i);
    imagesc(squeeze(explained_var_mean(:,:,i)));
    axis image; axis off;
    colorbar;
end
%%
explained_var_mean2 = mean(explained_var_all,3);
figure;
imagesc(squeeze(explained_var_mean2));
axis image; axis off;
caxis([0.8,1]);
colorbar;