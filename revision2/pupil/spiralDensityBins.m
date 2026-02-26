function count_sample = spiralDensityBins(T,data_folder,spirals_sort)
%%
pixSize = 0.01;                                                            % mm/pix
pixArea = pixSize^2;
ssp_index = get_ssp_index(data_folder);
hist_bin = 40;
mouseN = size(T,1);
count_sample = zeros(mouseN,size(spirals_sort,2));
for kk  = 1:mouseN
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %%
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    t = readNPY(fullfile(session_root,'svdTemporalComponents_corr.timestamps.npy'));
    %%
    for m = 1:size(spirals_sort,2)
        clear spirals_temp unique_spirals unique_spirals_unit
        spirals_temp = spirals_sort{kk,m};
        [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
        [lia,locb] = ismember(unique_spirals(:,1:2),ssp_index,'rows');
        unique_spirals = unique_spirals(lia,:);
        unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit = unique_spirals_unit./numel(t)*35*size(spirals_sort,2);                % spirals/(mm^2*s)
        count_sample(kk,m) = max(unique_spirals_unit(:));
    end
end
