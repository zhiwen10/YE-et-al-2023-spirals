function getScatterDataVsFftn(T,data_folder,save_folder,freq)
%%
freq_name = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
local_data_folder = fullfile(data_folder,...
    'spirals\spirals_fftn');
control_data_folder = fullfile(local_data_folder,...
    freq_name,'control_stats');
fftn_data_folder = fullfile(local_data_folder,...
    freq_name,'fftn_stats');
%%
ssp_index = get_ssp_index(data_folder);
%%
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
hist_bin = 40;
count_sample_control = zeros(15,10);
count_sample_permute = zeros(15,10);
for kk = 1:size(T,1)
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    control = load(fullfile(control_data_folder,[fname '_density.mat']));
    permute = load(fullfile(fftn_data_folder,[fname '_density.mat']));
    count = 1;
    for radius = 10:10:100
        clear unique_spirals_unit_control unique_spirals_unit_permute 
        spirals_control = control.spiral_density{count};
        spirals_permute = permute.spiral_density{count};
        [lia,locb] = ismember(spirals_control(:,1:2),ssp_index,'rows');
        spirals_control = spirals_control(lia,:);
        [lia1,locb1] = ismember(spirals_permute(:,1:2),ssp_index,'rows');
        spirals_permute = spirals_permute(lia1,:);
        unique_spirals_unit_control = spirals_control(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_control = unique_spirals_unit_control./control.frame_all(1)*35; % spirals/(mm^2*s)
        unique_spirals_unit_permute = spirals_permute(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
        unique_spirals_unit_permute = unique_spirals_unit_permute./permute.frame_all(1)*35; % spirals/(mm^2*s)
        % interp histgram counts 
        count_sample_control(kk,count) = max(unique_spirals_unit_control(:));
        count_sample_permute(kk,count) = max(unique_spirals_unit_permute(:));
        count = count+1;
    end
end
%%
%% percentage change
for radius = 1:10
    control_all = squeeze(count_sample_control(:,radius));
    permute_all = squeeze(count_sample_permute(:,radius));
    percentage_change(:,radius) = (control_all-permute_all)./permute_all;
end
%%
save(fullfile(save_folder,'spirals_fft_radius.mat'),...
    'count_sample_control','count_sample_permute','percentage_change');