function getFFTNSpiralsStats(T,freq,label,data_folder,save_folder)
 %%
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
spirals_folder = fullfile(data_folder,'spirals','spirals_fftn',...
    freq_folder,label);
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));
scale = 1;
BW = logical(projectedAtlas1);
[row,col] = find(BW);
brain_index = [col,row];
%%
frame_all = 0;
hist_bin = 40;
for kk = 1:size(T,1)
    %% session info
    clear spiralsT
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load data
    fname = [mn '_' tdb '_' num2str(en)];   
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));   % load atlas transformation matrix tform
    %%    
    load(fullfile(spirals_folder,[fname '.mat']));
    [spiralsT(:,1),spiralsT(:,2)] = transformPointsForward(tform,pwAll1(:,1),pwAll1(:,2));
    pwAll1(:,1:2) = round(spiralsT); 
    [lia,locb] = ismember(pwAll1(:,1:2),brain_index,'rows');
    pwAll1 = pwAll1(lia,:);
    %%
    clear spiral_density
    count = 1;
    for radius = 10:10:100
        clear unique_spirals
        spirals_temp = [];
        spirals_temp = pwAll1(pwAll1(:,3)==radius,:);
        [unique_spirals,scolor,low_color_bound,high_color_bound] = ...
            density_color_plot(spirals_temp,hist_bin);
        spiral_density{count} = unique_spirals;
        frame_all(count) = frame_count;
        count = count+1;
    end
    save(fullfile(save_folder,freq_folder,[label '_stats'],...
        [fname '_density.mat']), ...
        'spiral_density','frame_all');
end
