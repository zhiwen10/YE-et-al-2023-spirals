function getFFTNSpiralsMap(T,freq,label,data_folder,save_folder)
%%
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
spirals_folder = fullfile(data_folder,'spirals','spirals_fftn',...
    freq_folder,label);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));
%%
spirals_all = [];
frame_all = 0;
for kk = 1:size(T,1)
    %% session info
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% read data
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(spirals_folder,[fname '.mat']));
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(...
        tform,pwAll1(:,1),pwAll1(:,2));
    pwAll1(:,1:2) = round(spiralsT); 
    spirals_all = [spirals_all;pwAll1];
    frame_all = frame_all+frame_count;
end
%%
spirals_all(:,1:2) = round(spirals_all(:,1:2));
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
[lia,locb] = ismember(spirals_all(:,1:2),brain_index,'rows');
spirals_all = spirals_all(lia,:);
hist_bin = 40;
for radius = 10:10:100
    clear unique_spirals
    spirals_temp = [];
    spirals_temp = spirals_all(spirals_all(:,3)==radius,:);
    [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
        density_color_plot(spirals_temp,hist_bin);
    fname = ['histogram_radius_' num2str(radius) '.mat'];
    save(fullfile(save_folder,freq_folder,[label '_map'],fname),...
        'unique_spirals');
end