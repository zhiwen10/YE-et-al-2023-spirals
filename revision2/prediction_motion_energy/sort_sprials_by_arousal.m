function hs13g = sort_sprials_by_arousal(T,data_folder,save_folder)
% for each sprials in the raw spiral list, find a (and only one) spiral 
% in the nearby 100 pixel distance in the predicted frame. 
% If exist, then save both raw and preidcted spiral as a row. 

% In the final saved cell structure (left and right): 
% each cell is a single session
% Within each cell, each row is a spiral that has a nearby spiral in
% prediction

% Col  1       2        3      4          5      
%    raw_x    raw_y  radius  direction  frameN

% Col  6       7        8      9          10      
%    pdt_x    pdt_y  radius  direction  frameN

% Col  11           12
%    match   frame 2-8Hz amp
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);                                             % atlas brain boundry binary mask 
[row,col] = find(BW);
brain_index = [col,row];  
%%
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
T = T(logical(T.facevideo),:);
T([14,19],:) = [];
%%
n_bins = 18; 
low_freq_band = [0.05, 0.5];                                               % Phase-providing frequency (slow oscillation)
spirals_sort = phaseSpiralHistogram2(T,data_folder,n_bins,brain_index,low_freq_band);
count_sample = spiralDensityBins2(T,data_folder,spirals_sort);
%%
color1 = cbrewer2('seq','Greys',size(count_sample,1)+3);
color2 = '-b';
phase_bins = linspace(-pi, pi, n_bins+1);
phase_centers = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;
%%
color2 = {'g','r','c','m'};
ratio_all = {};
ratio_perm = {};
hs13g = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
count1 = 1;
for kk = [1,2,4] % 'THAL', 'STR'; 'MB'
    %%
    subplot(1,3,count1);
    clear indx T2
    indx = contains(T.Area,area(1,kk));
    T2 = T(indx,:);
    count_sample_temp = count_sample(indx,:);
    %%
    for session = 1:size(count_sample_temp,1)
        spiral = squeeze(count_sample_temp(session,:));     
        % plot(phase_centers,spiral,'color',color1(session+3,:));
        plot(phase_centers,spiral,'color',[0.5,0.5,0.5]);
        hold on;
    end
    xticks([-pi:pi/2:pi]);
    xticklabels({'-pi','-1/2*pi','0','1/2*pi','pi'});
    xlim([-pi,pi]);
    ylim([0,6]);
    yticks([0:2:6]);
    count_mean = mean(count_sample_temp,1);
    count_sem = std(count_sample_temp,[],1)./sqrt(size(count_sample_temp,1));
    hold on;
    shadedErrorBar(phase_centers, count_mean, count_sem, 'lineprops', color2{kk});
    xlabel('Phase of face motion');
    ylabel('Peak density (Centers/mm2*s)');
    count1 = count1+1;
end
%%
print(hs13g,fullfile(save_folder, 'FigS13g_spirals_arousal.pdf'),...
    '-dpdf', '-bestfit', '-painters');
