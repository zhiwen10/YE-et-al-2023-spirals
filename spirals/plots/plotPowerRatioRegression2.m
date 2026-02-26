function hs1f = plotPowerRatioRegression2(data_folder,save_folder)
%%
nameList2 = {'VISp','RSP','SSp','MOs'};
% load power ratio
load(fullfile(data_folder,'spirals','spirals_power_spectrum2','power_ratio_all_sessions.mat'));
% power ratio vs sprial density regression
load(fullfile(data_folder,'spirals','spirals_density',...
    'spiralDensityLinePerSession.mat'))
count_sample(count_sample<0) = 0;
max_density = max(count_sample,[],1);
%%
hs1f = figure('Renderer', 'painters', 'Position', [100 100 300 300]);
ax1 = subplot(1,1,1);
mdla = fitlm(max_density,power_ratio(3,:));
r2a = mdla.Rsquared.Adjusted;
plot(mdla);
ylim([0,50])
xlabel('Spiral density (spirals/mm2*s)');
ylabel('2-8 Hz power ratio (%)');
xlabel('');
xticks([0:5]);      
xticklabels(num2cell([0:5]));
legend('off');
print(hs1f, fullfile(save_folder,'FigS1f_spiral_density_vs_power'),'-dpdf', '-bestfit', '-painters');
