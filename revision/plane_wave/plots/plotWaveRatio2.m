function hs8m = plotWaveRatio2(data_folder,save_folder)
%% load plane wave index
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all2.mat'));
vxy_all = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
vxy_all2 = sum(vxy_all,3)./4;
angle_all = angle(vxy_all2);
amp_all = abs(vxy_all2);
for i = 1:15
    amp_temp = amp_all(i,:);
    ratio(i,1) = sum(amp_temp>0.4)./numel(amp_temp);
end
%% load spiral wave peak density
load(fullfile(data_folder,'spirals','spirals_density',...
    'spiralDensityLinePerSession.mat'))
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
spiral_peak_ratio = max_density./35;
%% load plane wave index new sessions
load(fullfile(data_folder,'revision','plane_wave','flow_mirror_all_newsession2.mat'));
vxy_all_new = cat(3,vxy_MO_left,vxy_MO_right,vxy_SSp_left,vxy_SSp_right);
vxy_all2_new = sum(vxy_all_new,3)./4;
angle_all = angle(vxy_all2_new);
amp_all_new = abs(vxy_all2_new);
for i = 1:4
    amp_temp = amp_all_new(i,:);
    ratio_new(i,1) = sum(amp_temp>0.4)./numel(amp_temp);
end
%% load spiral wave peak density new sessions
load(fullfile(data_folder,'revision','plane_wave',...
    'spiralDensityLinePerSession_new.mat'))
count_sample(count_sample<0) = 0;
max_density_new = max(count_sample,[],1);
spiral_peak_ratio_new = max_density_new./35;
%%
mean_g7_spiral = mean(spiral_peak_ratio(1:11));
std_g7_spiral = std(spiral_peak_ratio(1:11));
sem_g7_spiral = std_g7_spiral./sqrt(11);

mean_g8_spiral = mean(spiral_peak_ratio(12:15));
std_g8_spiral = std(spiral_peak_ratio(12:15));
sem_g8_spiral = std_g8_spiral./sqrt(4);

mean_g7_plane = mean(ratio(1:11));
std_g7_plane = std(ratio(1:11));
sem_g7_plane = std_g7_plane./sqrt(11);

mean_g8_plane = mean(ratio(12:15));
std_g8_plane = std(ratio(12:15));
sem_g8_plane = std_g8_plane./sqrt(4);
%%
hs8m = figure('Renderer', 'painters', 'Position', [50 50 250 250]);
ax3 = subplot(1,1,1); 
scatter(ratio(1:11),spiral_peak_ratio(1:11),[],'k');
hold on;
scatter(ratio(12:15),spiral_peak_ratio(12:15),[],'r');

hold on;
errorbar(mean_g7_plane,mean_g7_spiral,sem_g7_spiral,sem_g7_spiral,sem_g7_plane,sem_g7_plane,...
'marker','None','lineWidth',1.5,'color','k');
hold on;
errorbar(mean_g8_plane,mean_g8_spiral,sem_g8_spiral,sem_g8_spiral,sem_g8_plane,sem_g8_plane,...
'marker','None','lineWidth',1.5,'color','r');
hold on;
plot([0, 0.4],[0, 0.4],'k--');

xlim([0,0.4]);
ylim([0,0.4]);
xlabel('Plane wave ratio');
ylabel('Spiral wave ratio');
xticks([0:0.1:0.4]);
xticklabels({'0','10%','20%','30%','40%'});
yticks([0:0.1:0.4]);
yticklabels({'0','10%','20%','30%','40%'});

%%
[hh,pp] = ttest(ratio(12:15),spiral_peak_ratio(12:15)');
%%
[hh1,plane_p1] = ttest2(ratio(1:11),ratio(12:15));
[hh1,spiral_p1] = ttest2(spiral_peak_ratio(1:11),spiral_peak_ratio(12:15));
%%
print(hs8m, fullfile(save_folder,'FigS8m_wave_ratio.pdf'),...
    '-dpdf', '-bestfit', '-painters');
